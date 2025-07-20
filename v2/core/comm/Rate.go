package comm

import (
	"strconv"
	"strings"
)

const (
	ALLOC_STEPS     = 6
	LOG_MAX_PSEUDO  = 6
	MAX_FINE_BITS   = 6
	FINE_OFFSET     = 21
	LOG2_FRAC_TABLE = "0,8,13,16,19,21,23,24,26,27,28,29,30,31,32,32,33,34,34,35,36,36,37,37"
)

// 获取脉冲数
func GetPulses(i int) int {
	if i < 8 {
		return i
	}
	return (8 + (i & 7)) << ((i >> 3) - 1)
}

// 比特数转脉冲数
func Bits2Pulses(m *CeltMode, band int, LM int, bits int) int {
	cache := m.cache.bits
	cachePtr := m.cache.index[LM*m.nbEBands+band]

	lo := 0
	hi := int(cache[cachePtr])
	bits--

	for i := 0; i < LOG_MAX_PSEUDO; i++ {
		mid := (lo + hi + 1) >> 1
		if int(cache[cachePtr+mid]) >= bits {
			hi = mid
		} else {
			lo = mid
		}
	}

	loVal := 0
	if lo > 0 {
		loVal = int(cache[cachePtr+lo])
	}
	hiVal := int(cache[cachePtr+hi])

	if bits-loVal <= hiVal-bits {
		return lo
	}
	return hi
}

// 脉冲数转比特数
func Pulses2Bits(m *CeltMode, band int, LM int, pulses int) int {
	LM++
	if pulses == 0 {
		return 0
	}
	return int(m.cache.bits[m.cache.index[LM*m.nbEBands+band]+pulses]) + 1
}

// 插值比特到脉冲
func InterpBits2Pulses(m *CeltMode, start int, end int, skipStart int,
	bits1 []int, bits2 []int, thresh []int, cap []int, total int, balance *int,
	skipRsv int, intensity *int, intensityRsv int, dualStereo *int, dualStereoRsv int, bits []int,
	ebits []int, finePriority []int, C int, LM int, ec *EntropyCoder, encode int, prev int, signalBandwidth int) int {

	allocFloor := C << BITRES
	stereo := 0
	if C > 1 {
		stereo = 1
	}

	logM := LM << BITRES
	lo := 0
	hi := 1 << ALLOC_STEPS

	// 二分查找
	for i := 0; i < ALLOC_STEPS; i++ {
		mid := (lo + hi) >> 1
		psum := 0
		done := 0
		for j := end - 1; j >= start; j-- {
			tmp := bits1[j] + (mid * bits2[j] >> ALLOC_STEPS)
			if tmp >= thresh[j] || done != 0 {
				done = 1
				psum += min(tmp, cap[j])
			} else if tmp >= allocFloor {
				psum += allocFloor
			}
		}
		if psum > total {
			hi = mid
		} else {
			lo = mid
		}
	}

	// 计算比特分配
	psum := 0
	done := 0
	for j := end - 1; j >= start; j-- {
		tmp := bits1[j] + (lo * bits2[j] >> ALLOC_STEPS)
		if tmp < thresh[j] && done == 0 {
			if tmp >= allocFloor {
				tmp = allocFloor
			} else {
				tmp = 0
			}
		} else {
			done = 1
		}
		tmp = min(tmp, cap[j])
		bits[j] = tmp
		psum += tmp
	}

	// 决定跳过哪些频带
	codedBands := end
	for {
		j := codedBands - 1
		if j <= skipStart {
			total += skipRsv
			break
		}

		left := total - psum
		percoeff := left / (m.eBands[codedBands] - m.eBands[start])
		rem := max(left-(m.eBands[j]-m.eBands[start])*percoeff, 0)
		bandWidth := m.eBands[codedBands] - m.eBands[j]
		bandBits := bits[j] + percoeff*bandWidth + rem

		allocFloorPlus := allocFloor + (1 << BITRES)
		if bandBits >= max(thresh[j], allocFloorPlus) {
			if encode != 0 {
				threshold := (min(prev, 7) * bandWidth << LM << BITRES) >> 4
				if codedBands <= start+2 || (bandBits > threshold && j <= signalBandwidth) {
					ec.EncBitLogp(1, 1)
					break
				}
				ec.EncBitLogp(0, 1)
			} else if ec.DecBitLogp(1) != 0 {
				break
			}
			psum += 1 << BITRES
			bandBits -= 1 << BITRES
		}

		psum -= bits[j] + intensityRsv
		if intensityRsv > 0 {
			intensityRsv = log2FracTable[j-start]
		}
		psum += intensityRsv

		if bandBits >= allocFloor {
			bits[j] = allocFloor
			psum += allocFloor
		} else {
			bits[j] = 0
		}
		codedBands--
	}

	// 编码强度和双立体声参数
	if intensityRsv > 0 {
		if encode != 0 {
			*intensity = min(*intensity, codedBands)
			ec.EncUint(uint32(*intensity-start), uint32(codedBands+1-start))
		} else {
			*intensity = start + int(ec.DecUint(uint32(codedBands+1-start)))
		}
	} else {
		*intensity = 0
	}

	if *intensity <= start {
		total += dualStereoRsv
		dualStereoRsv = 0
	}
	if dualStereoRsv > 0 {
		if encode != 0 {
			ec.EncBitLogp(*dualStereo, 1)
		} else {
			*dualStereo = ec.DecBitLogp(1)
		}
	} else {
		*dualStereo = 0
	}

	// 分配剩余比特
	left := total - psum
	percoeff := left / (m.eBands[codedBands] - m.eBands[start])
	left -= (m.eBands[codedBands] - m.eBands[start]) * percoeff
	for j := start; j < codedBands; j++ {
		bits[j] += percoeff * (m.eBands[j+1] - m.eBands[j])
	}
	for j := start; j < codedBands; j++ {
		tmp := min(left, m.eBands[j+1]-m.eBands[j])
		bits[j] += tmp
		left -= tmp
	}

	// 精细比特分配
	currBalance := 0
	for j := start; j < codedBands; j++ {
		N0 := m.eBands[j+1] - m.eBands[j]
		N := N0 << LM
		bit := bits[j] + currBalance

		if N > 1 {
			excess := max(bit-cap[j], 0)
			bits[j] = bit - excess

			den := C * N
			if C == 2 && N > 2 && *dualStereo == 0 && j < *intensity {
				den++
			}

			NClogN := den * (m.logN[j] + logM)
			offset := (NClogN >> 1) - den*FINE_OFFSET

			if N == 2 {
				offset += den << BITRES >> 2
			}

			if bits[j]+offset < den*2<<BITRES {
				offset += NClogN >> 2
			} else if bits[j]+offset < den*3<<BITRES {
				offset += NClogN >> 3
			}

			ebits[j] = max(0, (bits[j] + offset + (den << (BITRES - 1))))
			ebits[j] = ebits[j] / den >> BITRES

			if C*ebits[j] > (bits[j] >> BITRES) {
				ebits[j] = bits[j] >> stereo >> BITRES
			}

			ebits[j] = min(ebits[j], MAX_FINE_BITS)
			finePriority[j] = 0
			if ebits[j]*(den<<BITRES) >= bits[j]+offset {
				finePriority[j] = 1
			}

			bits[j] -= C * ebits[j] << BITRES
		} else {
			excess := max(bit-C<<BITRES, 0)
			bits[j] = bit - excess
			ebits[j] = 0
			finePriority[j] = 1
		}

		if excess > 0 {
			extraFine := min(excess>>(stereo+BITRES), MAX_FINE_BITS-ebits[j])
			ebits[j] += extraFine
			extraBits := extraFine * C << BITRES
			if extraBits >= excess-currBalance {
				finePriority[j] = 1
			}
			excess -= extraBits
		}
		currBalance = excess
	}
	*balance = currBalance

	// 跳过频带的处理
	for j := codedBands; j < end; j++ {
		ebits[j] = bits[j] >> stereo >> BITRES
		bits[j] = 0
		finePriority[j] = 0
		if ebits[j] < 1 {
			finePriority[j] = 1
		}
	}

	return codedBands
}

// 计算比特分配
func ComputeAllocation(m *CeltMode, start int, end int, offsets []int, cap []int, allocTrim int,
	intensity *int, dualStereo *int, total int, balance *int, pulses []int,
	ebits []int, finePriority []int, C int, LM int, ec *EntropyCoder, encode int, prev int, signalBandwidth int) int {

	total = max(total, 0)
	len := m.nbEBands
	skipStart := start
	skipRsv := 0
	if total >= 1<<BITRES {
		skipRsv = 1 << BITRES
	}
	total -= skipRsv

	intensityRsv := 0
	dualStereoRsv := 0
	if C == 2 {
		intensityRsv = log2FracTable[end-start]
		if intensityRsv > total {
			intensityRsv = 0
		} else {
			total -= intensityRsv
			if total >= 1<<BITRES {
				dualStereoRsv = 1 << BITRES
			}
			total -= dualStereoRsv
		}
	}

	bits1 := make([]int, len)
	bits2 := make([]int, len)
	thresh := make([]int, len)
	trimOffset := make([]int, len)

	for j := start; j < end; j++ {
		thresh[j] = max(C<<BITRES, (3*(m.eBands[j+1]-m.eBands[j])<<LM<<BITRES)>>4)
		trimOffset[j] = C * (m.eBands[j+1] - m.eBands[j]) * (allocTrim - 5 - LM) * (end - j - 1) * (1 << (LM + BITRES)) >> 6
		if (m.eBands[j+1]-m.eBands[j])<<LM == 1 {
			trimOffset[j] -= C << BITRES
		}
	}

	lo := 1
	hi := m.nbAllocVectors - 1
	for lo <= hi {
		mid := (lo + hi) >> 1
		psum := 0
		done := 0
		for j := end - 1; j >= start; j-- {
			bitsj := C * (m.eBands[j+1] - m.eBands[j]) * m.allocVectors[mid*len+j] << LM >> 2
			if bitsj > 0 {
				bitsj = max(0, bitsj+trimOffset[j])
			}
			bitsj += offsets[j]

			if bitsj >= thresh[j] || done != 0 {
				done = 1
				psum += min(bitsj, cap[j])
			} else if bitsj >= C<<BITRES {
				psum += C << BITRES
			}
		}
		if psum > total {
			hi = mid - 1
		} else {
			lo = mid + 1
		}
	}

	hi = lo
	lo = hi - 1
	for j := start; j < end; j++ {
		N := m.eBands[j+1] - m.eBands[j]
		bits1j := C * N * m.allocVectors[lo*len+j] << LM >> 2
		bits2j := cap[j]
		if hi < m.nbAllocVectors {
			bits2j = C * N * m.allocVectors[hi*len+j] << LM >> 2
		}

		if bits1j > 0 {
			bits1j = max(0, bits1j+trimOffset[j])
		}
		if bits2j > 0 {
			bits2j = max(0, bits2j+trimOffset[j])
		}
		if lo > 0 {
			bits1j += offsets[j]
		}
		bits2j += offsets[j]
		if offsets[j] > 0 {
			skipStart = j
		}
		bits2j = max(0, bits2j-bits1j)
		bits1[j] = bits1j
		bits2[j] = bits2j
	}

	codedBands := InterpBits2Pulses(m, start, end, skipStart, bits1, bits2, thresh, cap,
		total, balance, skipRsv, intensity, intensityRsv, dualStereo, dualStereoRsv, pulses,
		ebits, finePriority, C, LM, ec, encode, prev, signalBandwidth)

	return codedBands
}

// 初始化log2分数表
var log2FracTable []int

func init() {
	// 解析LOG2_FRAC_TABLE字符串
	strs := strings.Split(LOG2_FRAC_TABLE, ",")
	log2FracTable = make([]int, len(strs))
	for i, s := range strs {
		val, _ := strconv.Atoi(s)
		log2FracTable[i] = val
	}
}

// 辅助函数
func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}
