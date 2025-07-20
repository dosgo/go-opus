package opus

import (
	"errors"
)

var freqTableQ16 = []int32{
	12111, 9804, 8235, 7100, 6239, 5565, 5022, 4575, 4202,
	3885, 3612, 3375, 3167, 2984, 2820, 2674, 2542, 2422,
	2313, 2214, 2123, 2038, 1961, 1889, 1822, 1760, 1702,
}

const (
	WIN_TYPE_FIRST_HALF  = 1 // 0到π/2的正弦窗
	WIN_TYPE_SECOND_HALF = 2 // π/2到π的正弦窗
)

// ApplySineWindow 对信号应用正弦窗函数
func ApplySineWindow(pxWin, px []int16, winType, length int) error {
	// 参数验证
	if winType != WIN_TYPE_FIRST_HALF && winType != WIN_TYPE_SECOND_HALF {
		return errors.New("无效的窗口类型")
	}
	if length < 16 || length > 120 || (length&3) != 0 {
		return errors.New("长度必须是16-120之间且是4的倍数")
	}
	if len(px) < length || len(pxWin) < length {
		return errors.New("输入/输出缓冲区太小")
	}

	// 获取对应的频率值 (Q16格式)
	k := (length >> 2) - 4
	if k < 0 || k > 26 {
		return errors.New("计算窗口索引失败")
	}
	fQ16 := freqTableQ16[k]

	// 计算余弦近似因子 (Q16格式)
	cQ16 := SMULWB(fQ16, -fQ16)

	// 初始化状态
	var S0_Q16, S1_Q16 int32
	if winType == WIN_TYPE_FIRST_HALF {
		// 从0开始
		S0_Q16 = 0
		// 近似sin(f)
		S1_Q16 = fQ16 + int32(length>>3) // 右移3位相当于除以8
	} else {
		// 从1开始
		S0_Q16 = 1 << 16
		// 近似cos(f)
		S1_Q16 = (1 << 16) + (cQ16 >> 1) + int32(length>>4) // 右移4位相当于除以16
	}

	// 使用递归公式计算正弦窗:
	// sin(n*f) = 2 * cos(f) * sin((n-1)*f) - sin((n-2)*f)
	for k := 0; k < length; k += 4 {
		// 第一个样本: (S0 + S1) / 2 * px[k]
		sum := (S0_Q16 + S1_Q16) >> 1
		pxWin[k] = SMULWB(sum, int32(px[k]))

		// 第二个样本: S1 * px[k+1]
		pxWin[k+1] = SMULWB(S1_Q16, int32(px[k+1]))

		// 更新状态 S0 = 2*S1*cos(f) - S0 + 1
		S0_Q16 = SMULWB(S1_Q16, cQ16) + (S1_Q16 << 1) - S0_Q16 + 1
		S0_Q16 = minInt32(S0_Q16, 1<<16)

		// 第三个样本: (S0 + S1) / 2 * px[k+2]
		sum = (S0_Q16 + S1_Q16) >> 1
		pxWin[k+2] = SMULWB(sum, int32(px[k+2]))

		// 第四个样本: S0 * px[k+3]
		pxWin[k+3] = SMULWB(S0_Q16, int32(px[k+3]))

		// 更新状态 S1 = 2*S0*cos(f) - S1
		S1_Q16 = SMULWB(S0_Q16, cQ16) + (S0_Q16 << 1) - S1_Q16
		S1_Q16 = minInt32(S1_Q16, 1<<16)
	}

	return nil
}

// SMULWB 乘加运算 (a * b) >> 16 (带饱和)
func SMULWB(a, b int32) int16 {
	result := int64(a) * int64(b) >> 16
	if result > 32767 {
		return 32767
	}
	if result < -32768 {
		return -32768
	}
	return int16(result)
}

func minInt32(a, b int32) int32 {
	if a < b {
		return a
	}
	return b
}
