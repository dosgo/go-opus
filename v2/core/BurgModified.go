package opus

import (
	"github.com/project/silk/CeltPitchXCorr"
	"github.com/project/silk/Inlines"
	"github.com/project/silk/SilkConstants"
	"github.com/project/silk/TuningParameters"
)

const (
	MAX_FRAME_SIZE           = 384 // (0.005 * 16000+16)*4
	BURG_QA                  = 25
	BURG_N_BITS_HEAD_ROOM    = 2
	BURG_MIN_RSHIFTS         = -16
	BURG_MAX_RSHIFTS         = 32 - BURG_QA
	BURG_FIND_LPC_COND_FAC_Q = 32 // Q32 format
)

func Silk_burg_modified(
	res_nrg *BoxedValueInt, // O 残差能量
	res_nrg_Q *BoxedValueInt, // O 残差能量Q值
	A_Q16 []int32, // O 预测系数 (长度 order)
	x []int16, // I 输入信号
	x_ptr int, // I 输入信号偏移
	minInvGain_Q30 int32, // I 最大预测增益的倒数
	subfr_length int, // I 子帧长度 (包含前导样本)
	nb_subfr int, // I 子帧数量
	D int, // I 阶数
) {
	// 局部变量声明
	var k, n, s, lz, rshifts, reached_max_gain int
	var C0, num, nrg, rc_Q31, invGain_Q30, Atmp_QA, Atmp1, tmp1, tmp2, x1, x2 int32
	var x_offset int
	C_first_row := make([]int32, SilkConstants.SILK_MAX_ORDER_LPC)
	C_last_row := make([]int32, SilkConstants.SILK_MAX_ORDER_LPC)
	Af_QA := make([]int32, SilkConstants.SILK_MAX_ORDER_LPC)
	CAf := make([]int32, SilkConstants.SILK_MAX_ORDER_LPC+1)
	CAb := make([]int32, SilkConstants.SILK_MAX_ORDER_LPC+1)
	xcorr := make([]int32, SilkConstants.SILK_MAX_ORDER_LPC)
	var C0_64 int64

	// 断言检查
	if subfr_length*nb_subfr > MAX_FRAME_SIZE {
		panic("subfr_length*nb_subfr > MAX_FRAME_SIZE")
	}

	// 计算自相关
	C0_64 = Inlines.Silk_inner_prod16_aligned_64(x, x_ptr, x, x_ptr, subfr_length*nb_subfr)
	lz = Inlines.Silk_CLZ64(C0_64)
	rshifts = 32 + 1 + BURG_N_BITS_HEAD_ROOM - lz
	if rshifts > BURG_MAX_RSHIFTS {
		rshifts = BURG_MAX_RSHIFTS
	}
	if rshifts < BURG_MIN_RSHIFTS {
		rshifts = BURG_MIN_RSHIFTS
	}

	if rshifts > 0 {
		C0 = int32(Inlines.Silk_RSHIFT64(C0_64, rshifts))
	} else {
		C0 = Inlines.Silk_LSHIFT32(int32(C0_64), -rshifts)
	}

	// 初始化相关矩阵
	CAb[0] = C0 + Inlines.Silk_SMMUL(TuningParameters.FIND_LPC_COND_FAC_Q, C0) + 1
	CAf[0] = CAb[0]

	// 清零第一行
	for i := range C_first_row {
		C_first_row[i] = 0
	}

	// 计算第一行相关
	if rshifts > 0 {
		for s = 0; s < nb_subfr; s++ {
			x_offset = x_ptr + s*subfr_length
			for n = 1; n < D+1; n++ {
				corr := Inlines.Silk_inner_prod16_aligned_64(x, x_offset, x, x_offset+n, subfr_length-n)
				C_first_row[n-1] += int32(Inlines.Silk_RSHIFT64(corr, rshifts))
			}
		}
	} else {
		for s = 0; s < nb_subfr; s++ {
			x_offset = x_ptr + s*subfr_length
			CeltPitchXCorr.Pitch_xcorr(x, x_offset, x, x_offset+1, xcorr, subfr_length-D, D)

			// 补充计算边缘样本
			for n := 1; n < D+1; n++ {
				var d int32
				for i := n + subfr_length - D; i < subfr_length; i++ {
					d = Inlines.MAC16_16(d, x[x_offset+i], x[x_offset+i-n])
				}
				xcorr[n-1] += d
			}

			// 合并结果
			for n = 1; n < D+1; n++ {
				C_first_row[n-1] += Inlines.Silk_LSHIFT32(xcorr[n-1], -rshifts)
			}
		}
	}

	// 复制到最后一行
	copy(C_last_row, C_first_row)

	// 重新初始化
	CAb[0] = C0 + Inlines.Silk_SMMUL(TuningParameters.FIND_LPC_COND_FAC_Q, C0) + 1
	CAf[0] = CAb[0]

	// 主循环
	invGain_Q30 = 1 << 30
	reached_max_gain = 0
	for n = 0; n < D; n++ {
		// 更新相关矩阵
		if rshifts > -2 {
			for s = 0; s < nb_subfr; s++ {
				x_offset = x_ptr + s*subfr_length
				x1 = -Inlines.Silk_LSHIFT32(int32(x[x_offset+n]), 16-rshifts)
				x2 = -Inlines.Silk_LSHIFT32(int32(x[x_offset+subfr_length-n-1]), 16-rshifts)
				tmp1 = Inlines.Silk_LSHIFT32(int32(x[x_offset+n]), BURG_QA-16)
				tmp2 = Inlines.Silk_LSHIFT32(int32(x[x_offset+subfr_length-n-1]), BURG_QA-16)

				for k = 0; k < n; k++ {
					C_first_row[k] = Inlines.Silk_SMLAWB(C_first_row[k], x1, x[x_offset+n-k-1])
					C_last_row[k] = Inlines.Silk_SMLAWB(C_last_row[k], x2, x[x_offset+subfr_length-n+k])
					Atmp_QA = Af_QA[k]
					tmp1 = Inlines.Silk_SMLAWB(tmp1, Atmp_QA, x[x_offset+n-k-1])
					tmp2 = Inlines.Silk_SMLAWB(tmp2, Atmp_QA, x[x_offset+subfr_length-n+k])
				}

				tmp1 = Inlines.Silk_LSHIFT32(-tmp1, 32-BURG_QA-rshifts)
				tmp2 = Inlines.Silk_LSHIFT32(-tmp2, 32-BURG_QA-rshifts)

				for k = 0; k <= n; k++ {
					CAf[k] = Inlines.Silk_SMLAWB(CAf[k], tmp1, x[x_offset+n-k])
					CAb[k] = Inlines.Silk_SMLAWB(CAb[k], tmp2, x[x_offset+subfr_length-n+k-1])
				}
			}
		} else {
			// 类似处理，但使用不同移位
			for s = 0; s < nb_subfr; s++ {
				x_offset = x_ptr + s*subfr_length
				x1 = -Inlines.Silk_LSHIFT32(int32(x[x_offset+n]), -rshifts)
				x2 = -Inlines.Silk_LSHIFT32(int32(x[x_offset+subfr_length-n-1]), -rshifts)
				tmp1 = Inlines.Silk_LSHIFT32(int32(x[x_offset+n]), 17)
				tmp2 = Inlines.Silk_LSHIFT32(int32(x[x_offset+subfr_length-n-1]), 17)

				for k = 0; k < n; k++ {
					C_first_row[k] = Inlines.Silk_MLA(C_first_row[k], x1, x[x_offset+n-k-1])
					C_last_row[k] = Inlines.Silk_MLA(C_last_row[k], x2, x[x_offset+subfr_length-n+k])
					Atmp1 = Inlines.Silk_RSHIFT_ROUND(Af_QA[k], BURG_QA-17)
					tmp1 = Inlines.Silk_MLA(tmp1, x[x_offset+n-k-1], Atmp1)
					tmp2 = Inlines.Silk_MLA(tmp2, x[x_offset+subfr_length-n+k], Atmp1)
				}

				tmp1 = -tmp1
				tmp2 = -tmp2

				for k = 0; k <= n; k++ {
					shifted_x1 := Inlines.Silk_LSHIFT32(int32(x[x_offset+n-k]), -rshifts-1)
					shifted_x2 := Inlines.Silk_LSHIFT32(int32(x[x_offset+subfr_length-n+k-1]), -rshifts-1)
					CAf[k] = Inlines.Silk_SMLAWW(CAf[k], tmp1, shifted_x1)
					CAb[k] = Inlines.Silk_SMLAWW(CAb[k], tmp2, shifted_x2)
				}
			}
		}

		// 计算反射系数分子分母
		tmp1 = C_first_row[n]
		tmp2 = C_last_row[n]
		num = 0
		nrg = Inlines.Silk_ADD32(CAb[0], CAf[0])

		for k = 0; k < n; k++ {
			Atmp_QA = Af_QA[k]
			lz = Inlines.Silk_CLZ32(Inlines.Silk_abs(Atmp_QA)) - 1
			if lz > 32-BURG_QA {
				lz = 32 - BURG_QA
			}
			Atmp1 = Inlines.Silk_LSHIFT32(Atmp_QA, lz)

			tmp1 = Inlines.Silk_ADD_LSHIFT32(tmp1, Inlines.Silk_SMMUL(C_last_row[n-k-1], Atmp1), 32-BURG_QA-lz)
			tmp2 = Inlines.Silk_ADD_LSHIFT32(tmp2, Inlines.Silk_SMMUL(C_first_row[n-k-1], Atmp1), 32-BURG_QA-lz)
			num = Inlines.Silk_ADD_LSHIFT32(num, Inlines.Silk_SMMUL(CAb[n-k], Atmp1), 32-BURG_QA-lz)
			nrg = Inlines.Silk_ADD_LSHIFT32(nrg, Inlines.Silk_SMMUL(
				Inlines.Silk_ADD32(CAb[k+1], CAf[k+1]), Atmp1), 32-BURG_QA-lz)
		}

		CAf[n+1] = tmp1
		CAb[n+1] = tmp2
		num = Inlines.Silk_ADD32(num, tmp2)
		num = Inlines.Silk_LSHIFT32(-num, 1)

		// 计算反射系数
		var rc_Q31 int32
		if Inlines.Silk_abs(num) < nrg {
			rc_Q31 = Inlines.Silk_DIV32_varQ(num, nrg, 31)
		} else {
			if num > 0 {
				rc_Q31 = 2147483647 // INT32_MAX
			} else {
				rc_Q31 = -2147483648 // INT32_MIN
			}
		}

		// 更新逆增益
		tmp1 = (1 << 30) - Inlines.Silk_SMMUL(rc_Q31, rc_Q31)
		tmp1 = Inlines.Silk_LSHIFT(Inlines.Silk_SMMUL(invGain_Q30, tmp1), 2)
		if tmp1 <= minInvGain_Q30 {
			// 达到最大增益限制
			tmp2 = (1 << 30) - Inlines.Silk_DIV32_varQ(minInvGain_Q30, invGain_Q30, 30)
			rc_Q31 = Inlines.Silk_SQRT_APPROX(tmp2)
			rc_Q31 = Inlines.Silk_RSHIFT32(rc_Q31+Inlines.Silk_DIV32(tmp2, rc_Q31), 1)
			rc_Q31 = Inlines.Silk_LSHIFT32(rc_Q31, 16)
			if num < 0 {
				rc_Q31 = -rc_Q31
			}
			invGain_Q30 = minInvGain_Q30
			reached_max_gain = 1
		} else {
			invGain_Q30 = tmp1
		}

		// 更新AR系数
		for k = 0; k < (n+1)>>1; k++ {
			tmp1 = Af_QA[k]
			tmp2 = Af_QA[n-k-1]
			Af_QA[k] = Inlines.Silk_ADD_LSHIFT32(tmp1, Inlines.Silk_SMMUL(tmp2, rc_Q31), 1)
			Af_QA[n-k-1] = Inlines.Silk_ADD_LSHIFT32(tmp2, Inlines.Silk_SMMUL(tmp1, rc_Q31), 1)
		}
		Af_QA[n] = Inlines.Silk_RSHIFT32(rc_Q31, 31-BURG_QA)

		// 检查是否达到最大增益
		if reached_max_gain != 0 {
			for k = n + 1; k < D; k++ {
				Af_QA[k] = 0
			}
			break
		}

		// 更新CAf和CAb
		for k = 0; k <= n+1; k++ {
			tmp1 = CAf[k]
			tmp2 = CAb[n-k+1]
			CAf[k] = Inlines.Silk_ADD_LSHIFT32(tmp1, Inlines.Silk_SMMUL(tmp2, rc_Q31), 1)
			CAb[n-k+1] = Inlines.Silk_ADD_LSHIFT32(tmp2, Inlines.Silk_SMMUL(tmp1, rc_Q31), 1)
		}
	}

	// 处理结果
	if reached_max_gain != 0 {
		// 缩放系数
		for k = 0; k < D; k++ {
			A_Q16[k] = -Inlines.Silk_RSHIFT_ROUND(Af_QA[k], BURG_QA-16)
		}

		// 减去前导样本能量
		if rshifts > 0 {
			for s = 0; s < nb_subfr; s++ {
				x_offset = x_ptr + s*subfr_length
				corr := Inlines.Silk_inner_prod16_aligned_64(x, x_offset, x, x_offset, D)
				C0 -= int32(Inlines.Silk_RSHIFT64(corr, rshifts))
			}
		} else {
			for s = 0; s < nb_subfr; s++ {
				x_offset = x_ptr + s*subfr_length
				corr := Inlines.Silk_inner_prod_self(x, x_offset, D)
				C0 -= Inlines.Silk_LSHIFT32(corr, -rshifts)
			}
		}

		// 近似残差能量
		res_nrg.Val = int(Inlines.Silk_LSHIFT(Inlines.Silk_SMMUL(invGain_Q30, C0), 2))
		res_nrg_Q.Val = -rshifts
	} else {
		// 计算残差能量
		nrg = CAf[0]
		tmp1 = 1 << 16
		for k = 0; k < D; k++ {
			Atmp1 := Inlines.Silk_RSHIFT_ROUND(Af_QA[k], BURG_QA-16)
			nrg = Inlines.Silk_SMLAWW(nrg, CAf[k+1], Atmp1)
			tmp1 = Inlines.Silk_SMLAWW(tmp1, Atmp1, Atmp1)
			A_Q16[k] = -Atmp1
		}
		res_nrg.Val = int(Inlines.Silk_SMLAWW(nrg, Inlines.Silk_SMMUL(TuningParameters.FIND_LPC_COND_FAC_Q, C0), -tmp1))
		res_nrg_Q.Val = -rshifts
	}
}
