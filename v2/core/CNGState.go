package opus

const (
	MAX_FRAME_LENGTH = 576
	MAX_LPC_ORDER    = 16
)

// CNGState 舒适噪声生成(Comfort Noise Generation)状态结构
type CNGState struct {
	CNG_exc_buf_Q14   [MAX_FRAME_LENGTH]int32 // 随机激励缓冲区 (Q14格式)
	CNG_smth_NLSF_Q15 [MAX_LPC_ORDER]int16    // 平滑后的NLSF系数 (Q15格式)
	CNG_synth_state   [MAX_LPC_ORDER]int32    // 合成滤波器状态
	CNG_smth_Gain_Q32 int32                   // 平滑后的增益值 (Q16存储，但计算用Q32)
	rand_seed         int32                   // 随机数生成器种子
	Fs_kHz            int                     // 采样率 (kHz)
}

// Reset 重置CNG状态到初始值
func (c *CNGState) Reset(lpcOrder int) {
	// 清零激励缓冲区
	for i := range c.CNG_exc_buf_Q14 {
		c.CNG_exc_buf_Q14[i] = 0
	}

	// 初始化平滑NLSF为线性分布
	NLSF_step_Q15 := int32(32767) / int32(lpcOrder+1) // MAX_INT16/(LPC_order+1)
	NLSF_acc_Q15 := int32(0)

	for i := 0; i < lpcOrder; i++ {
		NLSF_acc_Q15 += NLSF_step_Q15
		c.CNG_smth_NLSF_Q15[i] = int16(NLSF_acc_Q15)
	}

	// 清零合成滤波器状态
	for i := range c.CNG_synth_state {
		c.CNG_synth_state[i] = 0
	}

	// 重置其他状态
	c.CNG_smth_Gain_Q32 = 0
	c.rand_seed = 3176576 // 原始Java代码中的初始种子
	c.Fs_kHz = 0
}
