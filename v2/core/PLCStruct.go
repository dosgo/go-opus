package opus

const (
	LTP_ORDER     = 5
	MAX_LPC_ORDER = 16
)

// PLCState 包丢失隐藏(Packet Loss Concealment)状态结构
type PLCState struct {
	pitchL_Q8         int32                // 用于语音隐藏的基音周期 (Q8格式)
	LTPCoef_Q14       [LTP_ORDER]int16     // 语音隐藏的LTP系数 (Q14格式)
	prevLPC_Q12       [MAX_LPC_ORDER]int16 // 前一帧的LPC系数 (Q12格式)
	last_frame_lost   int                  // 前一帧是否丢失标志
	rand_seed         int32                // 非语音随机信号生成种子
	randScale_Q14     int16                // 非语音随机信号缩放因子 (Q14格式)
	conc_energy       int32                // 能量水平（用于PLC）
	conc_energy_shift int32                // 能量水平的移位值
	prevLTP_scale_Q14 int16                // 前一帧的LTP缩放因子 (Q14格式)
	prevGain_Q16      [2]int32             // 前一帧的增益 (Q16格式)
	Fs_kHz            int                  // 采样率 (kHz)
	nb_subfr          int                  // 子帧数量
	subfr_length      int                  // 子帧长度
}

// Reset 重置PLC状态到初始值
func (s *PLCState) Reset() {
	s.pitchL_Q8 = 0
	for i := range s.LTPCoef_Q14 {
		s.LTPCoef_Q14[i] = 0
	}
	for i := range s.prevLPC_Q12 {
		s.prevLPC_Q12[i] = 0
	}
	s.last_frame_lost = 0
	s.rand_seed = 0
	s.randScale_Q14 = 0
	s.conc_energy = 0
	s.conc_energy_shift = 0
	s.prevLTP_scale_Q14 = 0
	s.prevGain_Q16[0] = 1 << 16 // Q16格式的1.0
	s.prevGain_Q16[1] = 1 << 16 // Q16格式的1.0
	s.Fs_kHz = 0
	s.nb_subfr = 0
	s.subfr_length = 0
}
