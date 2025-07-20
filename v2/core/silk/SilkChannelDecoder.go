package silk

const (
	MAX_FRAME_LENGTH         = 576
	MAX_SUB_FRAME_LENGTH     = 144
	MAX_LPC_ORDER            = 16
	MAX_FRAMES_PER_PACKET    = 3
	MAX_NB_SUBFR             = 4
	LTP_MEM_LENGTH_MS        = 20
	SUB_FRAME_LENGTH_MS      = 5
	MIN_LPC_ORDER            = 10
	SHELL_CODEC_FRAME_LENGTH = 16
	FLAG_DECODE_NORMAL       = 0
	FLAG_DECODE_LBRR         = 1
	TYPE_NO_VOICE_ACTIVITY   = 0
)

type BoxedValueInt struct{ Val int }

// 解码器状态结构
type DecoderState struct {
	prev_gain_Q16           int32
	exc_Q14                 [MAX_FRAME_LENGTH]int32
	sLPC_Q14_buf            [MAX_LPC_ORDER]int32
	outBuf                  [MAX_FRAME_LENGTH + 2*MAX_SUB_FRAME_LENGTH]int16
	lagPrev                 int
	LastGainIndex           int
	Fs_kHz                  int // 采样频率(kHz)
	Fs_API_hz               int // API采样频率(Hz)
	nb_subfr                int // 子帧数量
	frame_length            int // 帧长度(样本)
	subfr_length            int // 子帧长度(样本)
	ltp_mem_length          int // LTP内存长度
	LPC_order               int // LPC阶数
	prevNLSF_Q15            [MAX_LPC_ORDER]int16
	first_frame_after_reset int

	// 熵编码表指针
	pitch_lag_low_bits_iCDF []uint16
	pitch_contour_iCDF      []uint16

	// 多帧数据包处理
	nFramesDecoded   int
	nFramesPerPacket int

	// 熵编码状态
	ec_prevSignalType int
	ec_prevLagIndex   int16

	VAD_flags  [MAX_FRAMES_PER_PACKET]int
	LBRR_flag  int
	LBRR_flags [MAX_FRAMES_PER_PACKET]int

	resampler_state ResamplerState
	psNLSF_CB       *NLSFCB // NLSF码本指针

	indices SideInfoIndices // 量化索引

	sCNG CNGState // CNG状态

	// PLC相关状态
	lossCnt        int
	prevSignalType int
	sPLC           PLCState
}

// 重置解码器状态
func (s *DecoderState) Reset() {
	s.prev_gain_Q16 = 0
	for i := range s.exc_Q14 {
		s.exc_Q14[i] = 0
	}
	for i := range s.sLPC_Q14_buf {
		s.sLPC_Q14_buf[i] = 0
	}
	for i := range s.outBuf {
		s.outBuf[i] = 0
	}
	s.lagPrev = 0
	s.LastGainIndex = 0
	s.Fs_kHz = 0
	s.Fs_API_hz = 0
	s.nb_subfr = 0
	s.frame_length = 0
	s.subfr_length = 0
	s.ltp_mem_length = 0
	s.LPC_order = 0
	for i := range s.prevNLSF_Q15 {
		s.prevNLSF_Q15[i] = 0
	}
	s.first_frame_after_reset = 0
	s.nFramesDecoded = 0
	s.nFramesPerPacket = 0
	s.ec_prevSignalType = 0
	s.ec_prevLagIndex = 0
	for i := range s.VAD_flags {
		s.VAD_flags[i] = 0
	}
	s.LBRR_flag = 0
	for i := range s.LBRR_flags {
		s.LBRR_flags[i] = 0
	}
	s.resampler_state.Reset()
	s.psNLSF_CB = nil
	s.indices.Reset()
	s.sCNG.Reset(s.LPC_order)
	s.lossCnt = 0
	s.prevSignalType = 0
	s.sPLC.Reset()
}

// 初始化解码器
func (s *DecoderState) Init() int {
	s.Reset()
	s.first_frame_after_reset = 1
	s.prev_gain_Q16 = 65536 // 1 << 16
	s.sCNG.Reset(s.LPC_order)
	s.sPLC.Reset()
	return 0
}

// 设置采样率
func (s *DecoderState) SetSampleRate(fs_kHz, fs_API_Hz int) int {
	// 验证输入参数
	if fs_kHz != 8 && fs_kHz != 12 && fs_kHz != 16 {
		return -1
	}
	if s.nb_subfr != MAX_NB_SUBFR && s.nb_subfr != MAX_NB_SUBFR/2 {
		return -1
	}

	// 计算新的帧长
	subfr_len := SUB_FRAME_LENGTH_MS * fs_kHz
	frame_len := s.nb_subfr * subfr_len

	// 需要重采样器初始化的情况
	if s.Fs_kHz != fs_kHz || s.Fs_API_hz != fs_API_Hz {
		if err := ResamplerInit(&s.resampler_state, fs_kHz*1000, fs_API_Hz, 0); err != nil {
			return err
		}
		s.Fs_API_hz = fs_API_Hz
	}

	// 需要更新状态的情况
	if s.Fs_kHz != fs_kHz || frame_len != s.frame_length {
		switch {
		case fs_kHz == 8:
			if s.nb_subfr == MAX_NB_SUBFR {
				s.pitch_contour_iCDF = pitch_contour_NB_iCDF
			} else {
				s.pitch_contour_iCDF = pitch_contour_10_ms_NB_iCDF
			}
		case s.nb_subfr == MAX_NB_SUBFR:
			s.pitch_contour_iCDF = pitch_contour_iCDF
		default:
			s.pitch_contour_iCDF = pitch_contour_10_ms_iCDF
		}

		if s.Fs_kHz != fs_kHz {
			s.ltp_mem_length = LTP_MEM_LENGTH_MS * fs_kHz
			switch {
			case fs_kHz == 8 || fs_kHz == 12:
				s.LPC_order = MIN_LPC_ORDER
				s.psNLSF_CB = &NLSF_CB_NB_MB
			default:
				s.LPC_order = MAX_LPC_ORDER
				s.psNLSF_CB = &NLSF_CB_WB
			}

			switch fs_kHz {
			case 16:
				s.pitch_lag_low_bits_iCDF = uniform8_iCDF
			case 12:
				s.pitch_lag_low_bits_iCDF = uniform6_iCDF
			case 8:
				s.pitch_lag_low_bits_iCDF = uniform4_iCDF
			default:
				return -1 // 不支持的采样率
			}

			s.first_frame_after_reset = 1
			s.lagPrev = 100
			s.LastGainIndex = 10
			s.prevSignalType = TYPE_NO_VOICE_ACTIVITY
			for i := range s.outBuf {
				s.outBuf[i] = 0
			}
			for i := range s.sLPC_Q14_buf {
				s.sLPC_Q14_buf[i] = 0
			}
		}

		s.Fs_kHz = fs_kHz
		s.frame_length = frame_len
		s.subfr_length = subfr_len
	}

	// 验证设置
	if s.frame_length <= 0 || s.frame_length > MAX_FRAME_LENGTH {
		return -1
	}

	return 0
}

// 解码音频帧
func (s *DecoderState) DecodeFrame(
	rangeDec *RangeCoder, // 范围解码器
	out []int16, // 输出音频帧
	lostFlag int, // 丢包标志
	condCoding int, // 条件编码类型
) (n int, ret int) {
	decCtrl := DecoderControl{}
	L := s.frame_length
	decCtrl.LTP_scale_Q14 = 0

	// 检查帧长度有效性
	if L <= 0 || L > MAX_FRAME_LENGTH {
		return 0, -1
	}

	// 正常或LBRR解码
	if lostFlag == FLAG_DECODE_NORMAL ||
		(lostFlag == FLAG_DECODE_LBRR && s.LBRR_flags[s.nFramesDecoded] == 1) {

		// 解码量化索引
		if err := DecodeIndices(s, rangeDec, s.nFramesDecoded, lostFlag, condCoding); err != nil {
			return 0, err
		}

		// 准备脉冲缓冲区 (填充为SHELL_CODEC_FRAME_LENGTH的倍数)
		pulseLength := (L + SHELL_CODEC_FRAME_LENGTH - 1) & ^(SHELL_CODEC_FRAME_LENGTH - 1)
		pulses := make([]int16, pulseLength)

		// 解码激励脉冲
		if err := DecodePulses(rangeDec, pulses, s.indices.signalType,
			s.indices.quantOffsetType, s.frame_length); err != nil {
			return 0, err
		}

		// 解码核心参数
		if err := DecodeParameters(s, &decCtrl, condCoding); err != nil {
			return 0, err
		}

		// 核心解码处理
		if err := DecodeCore(s, &decCtrl, out, pulses); err != nil {
			return 0, err
		}

		// 更新PLC状态
		PLC_Update(s, &decCtrl, out, 0)

		// 更新帧状态
		s.lossCnt = 0
		s.prevSignalType = s.indices.signalType
		s.first_frame_after_reset = 0
	} else {
		// PLC处理
		PLC_Update(s, &decCtrl, out, 1)
	}

	// 更新输出缓冲区 (移动历史数据)
	mv_len := s.ltp_mem_length - s.frame_length
	if mv_len > 0 {
		copy(s.outBuf[:], s.outBuf[s.frame_length:])
	}
	copy(s.outBuf[mv_len:], out)

	// CNG处理
	CNG_Process(s, &decCtrl, out, L)

	// 平滑拼接帧
	PLC_GlueFrames(s, out, L)

	// 更新解码器状态
	s.lagPrev = decCtrl.pitchL[s.nb_subfr-1]
	s.nFramesDecoded++

	return L, 0
}
