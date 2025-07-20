package comm

import (
	"math"
	"math/bits"
)

// 常量定义
const (
	SYM_BITS   = 8
	SYM_MAX    = (1 << SYM_BITS) - 1
	CODE_BITS  = 32
	CODE_SHIFT = CODE_BITS - SYM_BITS - 1
	CODE_TOP   = 1 << (CODE_BITS - 1)
	CODE_BOT   = CODE_TOP >> SYM_BITS
	CODE_EXTRA = (CODE_BITS-2)%SYM_BITS + 1
	UNI_BITS   = 8
)

// ICDFContext 表示累积分布函数上下文
type ICDFContext struct {
	Total int
	Dist  []int
}

// RangeDecoder Opus 范围解码器
type RangeDecoder struct {
	bits       *UnpaddedBitReadBE
	revs       *ReverseBitReadLE
	Range      int
	Value      int
	Total      int
	SizeInBits int
}

// NewRangeDecoder 创建新的范围解码器
func NewRangeDecoder(buf []byte) *RangeDecoder {
	bits := NewUnpaddedBitReadBE(buf)
	value := 127 - bits.GetBits32(7)

	rd := &RangeDecoder{
		bits:       bits,
		revs:       NewReverseBitReadLE(buf),
		Range:      128,
		Value:      int(value),
		Total:      SYM_BITS + 1,
		SizeInBits: len(buf) * 8,
	}

	rd.normalize()
	return rd
}

// normalize 规范化范围解码器状态
func (rd *RangeDecoder) normalize() {
	for rd.Range <= CODE_BOT {
		v := rd.bits.GetBits32(SYM_BITS)
		v = v ^ SYM_MAX
		rd.Value = ((rd.Value << SYM_BITS) | int(v)) & (CODE_TOP - 1)
		rd.Range <<= SYM_BITS
		rd.Total += SYM_BITS
	}
}

// update 更新解码器状态
func (rd *RangeDecoder) update(scale, low, high, total int) {
	s := scale * (total - high)
	rd.Value -= s
	if low != 0 {
		rd.Range = scale * (high - low)
	} else {
		rd.Range -= s
	}

	if rd.Range == 0 {
		panic("range became zero")
	}

	rd.normalize()
}

// getScaleSymbol 获取缩放因子和符号
func (rd *RangeDecoder) getScaleSymbol(total int) (int, int) {
	scale := rd.Range / total
	k := total - min(rd.Value/scale+1, total)
	return scale, k
}

// DecodeLogP 解码对数概率值
func (rd *RangeDecoder) DecodeLogP(logp int) bool {
	scale := rd.Range >> logp
	k := scale > rd.Value

	if k {
		rd.Range = scale
	} else {
		rd.Range -= scale
		rd.Value -= scale
	}

	rd.normalize()
	return k
}

// DecodeICDF 解码累积分布函数
func (rd *RangeDecoder) DecodeICDF(icdf *ICDFContext) int {
	total := icdf.Total
	dist := icdf.Dist
	scale, sym := rd.getScaleSymbol(total)

	k := 0
	for i, v := range dist {
		if v > sym {
			k = i
			break
		}
	}

	high := dist[k]
	low := 0
	if k > 0 {
		low = dist[k-1]
	}

	rd.update(scale, low, high, total)
	return k
}

// Tell 获取当前比特位置
func (rd *RangeDecoder) Tell() int {
	return rd.Total - CeltIlog2(rd.Range)
}

// TellFrac 获取分数比特位置
func (rd *RangeDecoder) TellFrac() int {
	lg := CeltIlog2(rd.Range)
	rq15 := rd.Range >> (lg - 16)

	for i := 0; i < 3; i++ {
		rq15 = (rq15 * rq15) >> 15
		lastbit := rq15 >> 16
		lg = lg*2 | int(lastbit)
		rq15 >>= lastbit
	}

	return rd.Total*8 - lg
}

// Len 获取总比特长度
func (rd *RangeDecoder) Len() int {
	return rd.SizeInBits
}

// Available 获取可用比特数
func (rd *RangeDecoder) Available() int {
	return rd.SizeInBits - rd.Tell()
}

// AvailableFrac 获取分数可用比特数
func (rd *RangeDecoder) AvailableFrac() int {
	return rd.SizeInBits*8 - rd.TellFrac()
}

// CeltOnly 接口定义 Celt 特定解码方法
type CeltOnly interface {
	RawBits(len int) int
	DecodeUniform(len int) int
	DecodeLaplace(symbol, decay int) int
	DecodeStep(k0 int) int
	DecodeTriangular(qn int) int
	ToEnd()
}

// RawBits 解码原始比特
func (rd *RangeDecoder) RawBits(len int) int {
	rd.Total += len
	return int(rd.revs.GetBits32(len))
}

// DecodeUniform 解码均匀分布值
func (rd *RangeDecoder) DecodeUniform(len int) int {
	bits := CeltIlog2(len - 1)
	total := len
	if bits > UNI_BITS {
		total = ((len - 1) >> (bits - UNI_BITS)) + 1
	}

	scale, k := rd.getScaleSymbol(total)
	rd.update(scale, k, k+1, total)

	if bits > UNI_BITS {
		return (k << (bits - UNI_BITS)) | rd.RawBits(bits-UNI_BITS)
	}
	return k
}

// DecodeLaplace 解码拉普拉斯分布值
func (rd *RangeDecoder) DecodeLaplace(symbol, decay int) int {
	scale := rd.Range >> 15
	center := rd.Value/scale + 1
	center = min(center, 1<<15)
	center = (1 << 15) - center

	value := 0
	low := 0

	if center >= symbol {
		value = 1
		low = symbol
		symbol = 1 + ((32768 - 32 - symbol) * (16384 - decay) >> 15)

		for symbol > 1 && center >= low+2*symbol {
			value++
			symbol *= 2
			low += symbol
			symbol = (((symbol - 2) * decay) >> 15) + 1
		}

		if symbol <= 1 {
			dist := (center - low) >> 1
			value += dist
			low += 2 * dist
		}

		if center < low+symbol {
			value = -value
		} else {
			low += symbol
		}
	}

	rd.update(scale, low, min(low+symbol, 32768), 32768)
	return value
}

// DecodeStep 解码步进值
func (rd *RangeDecoder) DecodeStep(k0 int) int {
	k1 := (k0 + 1) * 3
	total := k1 + k0
	scale := rd.Range / total
	symbol := rd.Value/scale + 1
	symbol = total - min(symbol, total)

	k := 0
	if symbol < k1 {
		k = symbol / 3
	} else {
		k = symbol - (k0+1)/2
	}

	if k <= k0 {
		rd.update(scale, 3*k, 3*(k+1), total)
	} else {
		rd.update(scale, 3*(k+1)+(k-1-k0), 3*(k0+1)+(k-k0), total)
	}

	return k
}

// DecodeTriangular 解码三角分布值
func (rd *RangeDecoder) DecodeTriangular(qn int) int {
	qn2 := qn >> 1
	total := (qn2 + 1) * (qn2 + 1)
	scale := rd.Range / total
	center := rd.Value/scale + 1
	center = total - min(center, total)

	k := 0
	low := 0
	symbol := 0

	if center < total>>1 {
		k = int((math.Sqrt(float64(8*center+1)) - 1) / 2)
		low = k * (k + 1) >> 1
		symbol = k + 1
	} else {
		k = (2*(qn+1) - int(math.Sqrt(float64(8*(total-center-1)+1)))>>1)
		low = total - ((qn + 1 - k) * (qn + 2 - k) >> 1)
		symbol = qn + 1 - k
	}

	rd.update(scale, low, low+symbol, total)
	return k
}

// ToEnd 移动到比特流末尾
func (rd *RangeDecoder) ToEnd() {
	rd.Total += rd.SizeInBits - rd.Tell()
}

// CeltIlog2 计算以2为底的对数
func CeltIlog2(x int) int {
	if x <= 0 {
		return 0
	}
	return 31 - bits.LeadingZeros32(uint32(x))
}

// min 返回两个整数中的最小值
func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

// 以下是辅助的比特读取器实现

// UnpaddedBitReadBE 大端序无填充比特读取器
type UnpaddedBitReadBE struct {
	buffer []byte
	index  int
}

// NewUnpaddedBitReadBE 创建新的大端序比特读取器
func NewUnpaddedBitReadBE(buf []byte) *UnpaddedBitReadBE {
	return &UnpaddedBitReadBE{
		buffer: buf,
		index:  0,
	}
}

// GetBits32 读取指定长度的比特
func (r *UnpaddedBitReadBE) GetBits32(bits int) uint32 {
	value := uint32(0)
	for i := 0; i < bits; i++ {
		if r.index < len(r.buffer) {
			value = (value << 1) | uint32((r.buffer[r.index]>>(7-i%8))&1)
			if (i+1)%8 == 0 {
				r.index++
			}
		}
	}
	return value
}

// ReverseBitReadLE 小端序反向比特读取器
type ReverseBitReadLE struct {
	buffer []byte
	index  int
}

// NewReverseBitReadLE 创建新的小端序反向比特读取器
func NewReverseBitReadLE(buf []byte) *ReverseBitReadLE {
	return &ReverseBitReadLE{
		buffer: buf,
		index:  len(buf) - 1,
	}
}

// GetBits32 读取指定长度的比特
func (r *ReverseBitReadLE) GetBits32(bits int) uint32 {
	value := uint32(0)
	for i := 0; i < bits; i++ {
		if r.index >= 0 {
			value = (value << 1) | uint32((r.buffer[r.index]>>(i%8))&1)
			if (i+1)%8 == 0 {
				r.index--
			}
		}
	}
	return value
}
