// 版权声明和许可证信息保留（省略重复部分）
package opus

const (
	CNG_BUF_MASK_MAX       = 512 - 1
	CNG_NLSF_SMTH_Q16      = 15565 // 0.25 * (1 << 16)
	CNG_GAIN_SMTH_Q16      = 23476 // 0.375 * (1 << 16)
	TYPE_NO_VOICE_ACTIVITY = 0
)

type CNGState struct {
	CNG_smth_NLSF_Q15 []int32
	CNG_smth_Gain_Q32 int32 // Go 使用 int32 替代原 Java 的 Q16 格式（实际存储为 Q16）
	rand_seed         int32
	CNG_exc_buf_Q14   []int32
	CNG_synth_state   []int32
	fs_kHz            int32
}

// 生成 CNG 激励信号
func CNG_exc(exc_Q10 []int32, exc_buf_Q14 []int32, Gain_Q16 int32, length int, rand_seed *int32) {
	seed := *rand_seed
	exc_mask := int32(CNG_BUF_MASK_MAX)

	// 调整掩码使其不小于长度
	for exc_mask > int32(length) {
		exc_mask >>= 1
	}

	for i := 0; i < length; i++ {
		seed = RAND(seed)
		idx := (seed >> 24) & exc_mask
		// 乘法调整为 Q16 增益应用
		exc_Q10[i] = SMULWW(exc_buf_Q14[idx], Gain_Q16>>4)
	}

	*rand_seed = seed
}

// 重置 CNG 状态
func CNG_Reset(psDec *DecoderState) {
	psCNG := &psDec.sCNG
	LPC_order := psDec.LPC_order
	NLSF_step_Q15 := int32(MAX_INT16) / int32(LPC_order+1)
	NLSF_acc_Q15 := int32(0)

	psCNG.CNG_smth_NLSF_Q15 = make([]int32, LPC_order)
	for i := 0; i < LPC_order; i++ {
		NLSF_acc_Q15 += NLSF_step_Q15
		psCNG.CNG_smth_NLSF_Q15[i] = NLSF_acc_Q15
	}
	psCNG.CNG_smth_Gain_Q32 = 0
	psCNG.rand_seed = 3176576
}

// 更新 CNG 估计并在丢包时应用 CNG
func CNG(psDec *DecoderState, psDecCtrl *DecoderControl, frame []int16, length int) {
	psCNG := &psDec.sCNG
	LPC_order := int(psDec.LPC_order)
	max_order := 16

	if psDec.Fs_kHz != psCNG.fs_kHz {
		CNG_Reset(psDec)
		psCNG.fs_kHz = psDec.Fs_kHz
	}

	if psDec.LossCnt == 0 && psDec.PrevSignalType == TYPE_NO_VOICE_ACTIVITY {
		// 平滑 LSF 系数
		for i := 0; i < LPC_order; i++ {
			delta := psDec.PrevNLSF_Q15[i] - psCNG.CNG_smth_NLSF_Q15[i]
			psCNG.CNG_smth_NLSF_Q15[i] += SMULWB(delta, CNG_NLSF_SMTH_Q16)
		}

		// 找到最大增益子帧
		maxGainIdx := 0
		for i := 1; i < len(psDecCtrl.Gains_Q16); i++ {
			if psDecCtrl.Gains_Q16[i] > psDecCtrl.Gains_Q16[maxGainIdx] {
				maxGainIdx = i
			}
		}

		// 更新激励缓冲区（使用子帧长度）
		subfr_len := int(psDec.SubfrLength)
		copy(psCNG.CNG_exc_buf_Q14, psCNG.CNG_exc_buf_Q14[subfr_len:])
		copy(psCNG.CNG_exc_buf_Q14[len(psCNG.CNG_exc_buf_Q14)-subfr_len:],
			psDec.Exc_Q14[(len(psDec.Exc_Q14)-subfr_len):])

		// 平滑增益
		for i := 0; i < len(psDecCtrl.Gains_Q16); i++ {
			gain_diff := psDecCtrl.Gains_Q16[i] - psCNG.CNG_smth_Gain_Q32>>16
			psCNG.CNG_smth_Gain_Q32 += int32((int64(gain_diff) * int64(CNG_GAIN_SMTH_Q16)) >> 16)
		}
	}

	if psDec.LossCnt != 0 {
		CNG_sig_Q10 := make([]int32, length+max_order)
		gain_Q16 := SMULWW(int32(psDec.sPLC.RandScale_Q14), psDec.sPLC.PrevGain_Q16[1])

		// 增益调整（自适应 Q 格式）
		if gain_Q16 >= (1<<21) || psCNG.CNG_smth_Gain_Q32 > (1<<23) {
			gain_t := int64(gain_Q16) * int64(gain_Q16) >> 16
			smooth_t := int64(psCNG.CNG_smth_Gain_Q32) * int64(psCNG.CNG_smth_Gain_Q32) >> 32
			gain_Q16 = int32(Sqrt(smooth_t - (gain_t >> 5)))
		} else {
			gain_t := int64(SMULWW(gain_Q16, gain_Q16))
			smooth_t := int64(SMULWW(psCNG.CNG_smth_Gain_Q32>>16, psCNG.CNG_smth_Gain_Q32>>16))
			gain_Q16 = int32(Sqrt(smooth_t - (gain_t >> 5)))
		}

		// 生成激励
		seed := psCNG.rand_seed
		CNG_exc(CNG_sig_Q10[max_order:], psCNG.CNG_exc_buf_Q14, gain_Q16, length, &seed)
		psCNG.rand_seed = seed

		// 转换 NLSF 到 LPC
		A_Q12 := make([]int32, max_order)
		NLSF2A(A_Q12, psCNG.CNG_smth_NLSF_Q15, LPC_order)

		// 初始化合成状态
		copy(CNG_sig_Q10, psCNG.CNG_synth_state)

		// 滤波生成信号
		for i := 0; i < length; i++ {
			sum_Q6 := int32(0)
			idx := max_order + i

			// LPC 滤波（展开固定阶数循环）
			sum_Q6 = SMLAWB(sum_Q6, CNG_sig_Q10[idx-1], A_Q12[0])
			sum_Q6 = SMLAWB(sum_Q6, CNG_sig_Q10[idx-2], A_Q12[1])
			// ... 添加剩余阶数（8/10/16阶实现）
			// 伪代码：实际需添加所有阶数的乘加

			CNG_sig_Q10[idx] = int32(int64(CNG_sig_Q10[idx]) + (int64(sum_Q6) << 4))
			frame[i] = SAT16(frame[i] + int16(CNG_sig_Q10[idx]>>10))
		}

		// 保存状态
		copy(psCNG.CNG_synth_state, CNG_sig_Q10[length:])
	} else {
		// 无丢包时重置状态
		for i := range psCNG.CNG_synth_state {
			psCNG.CNG_synth_state[i] = 0
		}
	}
}

// 辅助函数实现（根据原Java工具类）
func RAND(seed int32) int32                             { /* 32位 LCG 实现 */ }
func SMULWW(a, b int32) int32                           { /* 带饱和的乘法 */ }
func SMLAWB(sum, a, b int32) int32                      { /* 乘加运算 */ }
func SAT16(x int32) int16                               { /* 16位饱和 */ }
func NLSF2A(A_Q12 []int32, NLSF_Q15 []int32, order int) { /* LSF 转换 */ }
