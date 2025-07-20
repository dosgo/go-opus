package opus

import (
	"math"
)

const (
	QC = 10 // 扭曲自相关 Q 格式常数
	QS = 14 // 扭曲自相关 Q 格式常数
)

// Autocorr 计算标准自相关函数
func Autocorr(
	results []int32, // 输出: 自相关结果 (长度 correlationCount)
	scale *int, // 输出: 相关向量的缩放因子
	inputData []int16, // 输入: 待相关数据
	inputDataSize int, // 输入: 输入数据长度
	correlationCount int, // 输入: 要计算的相关点数
) {
	// 计算实际相关点数
	corrCount := minInt(inputDataSize, correlationCount)

	// 调用内部函数计算自相关
	*scale = celtAutocorr(inputData, results, corrCount-1, inputDataSize)
}

// celtAutocorr 自相关核心计算 (内部函数)
func celtAutocorr(
	x []int16, // 输入: [0...n-1] 样本
	ac []int32, // 输出: [0...lag-1] 自相关值
	lag int, // 输入: 最大滞后
	n int, // 输入: 样本数
) int {
	// 参数验证
	if n <= 0 {
		return 0
	}

	fastN := n - lag
	shift := 0
	xx := make([]int16, n)
	xptr := x

	// 计算初始能量并确定缩放因子
	ac0 := int32(1 + (n << 7))
	if n&1 != 0 {
		ac0 += MULT16_16(xptr[0], xptr[0]) >> 9
	}
	for i := n & 1; i < n; i += 2 {
		ac0 += MULT16_16(xptr[i], xptr[i]) >> 9
		ac0 += MULT16_16(xptr[i+1], xptr[i+1]) >> 9
	}

	// 计算缩放因子
	shift = celtIlog2(ac0) - 30 + 10
	shift = shift / 2

	// 应用缩放
	if shift > 0 {
		for i := 0; i < n; i++ {
			xx[i] = int16(PSHR32(int32(x[i]), shift))
		}
		xptr = xx
	} else {
		shift = 0
	}

	// 快速路径计算自相关
	acTmp := make([]int32, lag+1)
	PitchXCorr(xptr, xptr, acTmp, fastN, lag+1)

	// 慢速路径补充计算
	for k := 0; k <= lag; k++ {
		d := int32(0)
		for i := k + fastN; i < n; i++ {
			d = MAC16_16(d, xptr[i], xptr[i-k])
		}
		acTmp[k] += d
	}

	// 最终缩放调整
	shift *= 2
	if shift <= 0 {
		acTmp[0] += int32(1 << (-shift))
	}

	// 能量归一化
	if acTmp[0] < 268435456 { // 2^28
		shift2 := 29 - EC_ILOG(acTmp[0])
		for i := range acTmp {
			acTmp[i] <<= shift2
		}
		shift -= shift2
	} else if acTmp[0] >= 536870912 { // 2^29
		shift2 := 1
		if acTmp[0] >= 1073741824 { // 2^30
			shift2++
		}
		for i := range acTmp {
			acTmp[i] >>= shift2
		}
		shift += shift2
	}

	// 复制结果
	copy(ac, acTmp[:len(ac)])
	return shift
}

// WarpedAutocorrelation 计算扭曲频率轴的自相关
func WarpedAutocorrelation(
	corr []int32, // 输出: 结果 [order + 1]
	scale *int, // 输出: 相关向量的缩放因子
	input []int16, // 输入: 输入数据
	warpingQ16 int32, // 输入: 扭曲系数 (Q16格式)
	length int, // 输入: 输入长度
	order int, // 输入: 相关阶数 (偶数)
) {
	// 参数验证
	if order&1 != 0 {
		panic("阶数必须是偶数")
	}
	if 2*QS-QC < 0 {
		panic("QS/QC常数不兼容")
	}

	// 初始化状态
	stateQS := make([]int32, order+1)
	corrQC := make([]int64, order+1)

	// 处理每个样本
	for n := 0; n < length; n++ {
		tmp1QS := int32(input[n]) << QS

		// 处理全通滤波器段
		for i := 0; i < order; i += 2 {
			// 第一级全通输出
			tmp2QS := SMLAWB(stateQS[i], stateQS[i+1]-tmp1QS, warpingQ16)
			stateQS[i] = tmp1QS
			corrQC[i] += SMULL(tmp1QS, stateQS[0]) >> (2*QS - QC)

			// 第二级全通输出
			tmp1QS = SMLAWB(stateQS[i+1], stateQS[i+2]-tmp2QS, warpingQ16)
			stateQS[i+1] = tmp2QS
			corrQC[i+1] += SMULL(tmp2QS, stateQS[0]) >> (2*QS - QC)
		}

		// 更新最终状态
		stateQS[order] = tmp1QS
		corrQC[order] += SMULL(tmp1QS, stateQS[0]) >> (2*QS - QC)
	}

	// 计算缩放因子
	lsh := CLZ64(corrQC[0]) - 35
	lsh = clampInt(lsh, -12-QC, 30-QC)
	*scale = -(QC + lsh)

	// 应用缩放
	if lsh >= 0 {
		for i := range corr {
			corr[i] = int32(corrQC[i] << lsh)
		}
	} else {
		for i := range corr {
			corr[i] = int32(corrQC[i] >> (-lsh))
		}
	}
}

// ======== 定点运算辅助函数 ========

// MULT16_16 16位乘法 (a * b)
func MULT16_16(a, b int16) int32 {
	return int32(a) * int32(b)
}

// MAC16_16 乘加运算 (acc + a * b)
func MAC16_16(acc int32, a, b int16) int32 {
	return acc + int32(a)*int32(b)
}

// SMULL 64位乘法 (a * b)
func SMULL(a, b int32) int64 {
	return int64(a) * int64(b)
}

// SMLAWB 乘加运算 (a + (b * c) >> 16)
func SMLAWB(a, b, c int32) int32 {
	return a + (b*c)>>16
}

// PSHR32 带符号右移
func PSHR32(a int32, shift int) int32 {
	if shift > 0 {
		return a >> shift
	}
	return a << (-shift)
}

// celtIlog2 整数log2计算
func celtIlog2(x int32) int {
	return 31 - CLZ(uint32(x))
}

// EC_ILOG 有效位计算
func EC_ILOG(x int32) int {
	return 32 - CLZ(uint32(x))
}

// CLZ 计算前导零数量
func CLZ(x uint32) int {
	if x == 0 {
		return 32
	}
	return 31 - int(math.Log2(float64(x)))
}

// CLZ64 64位前导零计算
func CLZ64(x int64) int {
	if x == 0 {
		return 64
	}
	return 63 - int(math.Log2(float64(x)))
}

// clampInt 整数钳位
func clampInt(v, min, max int) int {
	if v < min {
		return min
	}
	if v > max {
		return max
	}
	return v
}

// minInt 整数最小值
func minInt(a, b int) int {
	if a < b {
		return a
	}
	return b
}
