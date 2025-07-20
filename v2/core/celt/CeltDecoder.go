package celt

import (
	"errors"
)

type CeltDecoder struct {
	mode           *CeltMode
	overlap        int
	channels       int
	streamChannels int
	downsample     int
	start          int
	end            int
	signalling     int

	// State cleared on reset
	rng                 int
	error               int
	lastPitchIndex      int
	lossCount           int
	postfilterPeriod    int
	postfilterPeriodOld int
	postfilterGain      int
	postfilterGainOld   int
	postfilterTapset    int
	postfilterTapsetOld int
	preemphMemD         [2]int

	// Dynamic buffers
	decodeMem      [][]int
	lpc            [][]int
	oldEBands      []int
	oldLogE        []int
	oldLogE2       []int
	backgroundLogE []int
}

func NewCeltDecoder() *CeltDecoder {
	return &CeltDecoder{
		preemphMemD: [2]int{0, 0},
	}
}

func (d *CeltDecoder) Reset() {
	d.mode = nil
	d.overlap = 0
	d.channels = 0
	d.streamChannels = 0
	d.downsample = 0
	d.start = 0
	d.end = 0
	d.signalling = 0
	d.PartialReset()
}

func (d *CeltDecoder) PartialReset() {
	d.rng = 0
	d.error = 0
	d.lastPitchIndex = 0
	d.lossCount = 0
	d.postfilterPeriod = 0
	d.postfilterPeriodOld = 0
	d.postfilterGain = 0
	d.postfilterGainOld = 0
	d.postfilterTapset = 0
	d.postfilterTapsetOld = 0
	d.preemphMemD = [2]int{0, 0}
	d.decodeMem = nil
	d.lpc = nil
	d.oldEBands = nil
	d.oldLogE = nil
	d.oldLogE2 = nil
	d.backgroundLogE = nil
}

func (d *CeltDecoder) ResetState() {
	d.PartialReset()

	if d.channels == 0 || d.mode == nil {
		return
	}

	d.decodeMem = make([][]int, d.channels)
	d.lpc = make([][]int, d.channels)
	for c := 0; c < d.channels; c++ {
		d.decodeMem[c] = make([]int, DECODE_BUFFER_SIZE+d.mode.overlap)
		d.lpc[c] = make([]int, LPC_ORDER)
	}

	nbEBands := d.mode.nbEBands
	d.oldEBands = make([]int, 2*nbEBands)
	d.oldLogE = make([]int, 2*nbEBands)
	d.oldLogE2 = make([]int, 2*nbEBands)
	d.backgroundLogE = make([]int, 2*nbEBands)

	for i := 0; i < 2*nbEBands; i++ {
		d.oldLogE[i] = -int(28.0 * float32(1<<DB_SHIFT))
		d.oldLogE2[i] = d.oldLogE[i]
	}
}

func (d *CeltDecoder) Init(samplingRate int, channels int) error {
	if err := d.OpusCustomDecoderInit(mode48000_960_120, channels); err != nil {
		return err
	}

	d.downsample = ResamplingFactor(samplingRate)
	if d.downsample == 0 {
		return errors.New("OPUS_BAD_ARG")
	}
	return nil
}

func (d *CeltDecoder) OpusCustomDecoderInit(mode *CeltMode, channels int) error {
	if channels < 1 || channels > 2 {
		return errors.New("OPUS_BAD_ARG")
	}

	d.Reset()
	d.mode = mode
	d.overlap = mode.overlap
	d.streamChannels = channels
	d.channels = channels
	d.downsample = 1
	d.start = 0
	d.end = mode.effEBands
	d.signalling = 1
	d.lossCount = 0
	d.ResetState()
	return nil
}

func (d *CeltDecoder) DecodeLost(N int, LM int) {
	// Implementation of PLC (Packet Loss Concealment)
	// This is a simplified version, actual implementation would require
	// pitch search, LPC analysis, and signal extrapolation
	for c := 0; c < d.channels; c++ {
		// Apply simple fade-out to avoid discontinuities
		for i := 0; i < N; i++ {
			fade := 1.0 - float64(i)/float64(N)
			d.decodeMem[c][DECODE_BUFFER_SIZE-N+i] = int(float64(d.decodeMem[c][DECODE_BUFFER_SIZE-N+i]) * fade)
		}
	}
	d.lossCount++
}

func (d *CeltDecoder) DecodeWithEC(data []byte, pcm []int16, frameSize int, ec *EntropyCoder, accum int) (int, error) {
	// 初始化变量
	var (
		c, i, N, LM, M, start, end, effEnd int
		spreadDecision, bits               int
		shortBlocks, isTransient           int
		intraEnerg, silence                int
		allocTrim, postfilterPitch         int
		postfilterGain, intensity          int
		dualStereo, totalBits, balance     int
		tell, dynallocLogp                 int
		postfilterTapset, antiCollapseRsv  int
		antiCollapseOn                     int
		C                                  = d.streamChannels
		CC                                 = d.channels
		mode                               = d.mode
	)

	// 计算帧大小
	frameSize *= d.downsample

	// 获取模式参数
	nbEBands := mode.nbEBands
	overlap := mode.overlap
	eBands := mode.eBands
	start = d.start
	end = d.end

	// 获取状态缓冲区
	oldBandE := d.oldEBands
	oldLogE := d.oldLogE
	oldLogE2 := d.oldLogE2
	backgroundLogE := d.backgroundLogE

	// 确定LM值
	for LM = 0; LM <= mode.maxLM; LM++ {
		if mode.shortMdctSize<<LM == frameSize {
			break
		}
	}
	if LM > mode.maxLM {
		return 0, errors.New("OPUS_BAD_ARG")
	}
	M = 1 << LM

	// 验证输入
	if len(data) > 1275 || pcm == nil {
		return 0, errors.New("OPUS_BAD_ARG")
	}

	// 计算N值
	N = M * mode.shortMdctSize

	// 准备输出缓冲区
	outSyn := make([][]int, CC)
	outSynPtrs := make([]int, CC)
	for c := 0; c < CC; c++ {
		outSyn[c] = d.decodeMem[c]
		outSynPtrs[c] = DECODE_BUFFER_SIZE - N
	}

	// 计算有效结束频带
	effEnd = end
	if effEnd > mode.effEBands {
		effEnd = mode.effEBands
	}

	// 处理丢包情况
	if len(data) <= 1 {
		d.DecodeLost(N, LM)
		Deemphasis(outSyn, outSynPtrs, pcm, N, CC, d.downsample, mode.preemph, d.preemphMemD[:], accum)
		return frameSize / d.downsample, nil
	}

	// 创建熵解码器（如果需要）
	var dec *EntropyCoder
	if ec == nil {
		dec = NewEntropyCoder(data)
	} else {
		dec = ec
	}

	// 单声道处理
	if C == 1 {
		for i := 0; i < nbEBands; i++ {
			if oldBandE[i] < oldBandE[nbEBands+i] {
				oldBandE[i] = oldBandE[nbEBands+i]
			}
		}
	}

	// 计算总比特数
	totalBits = len(data) * 8
	tell = dec.Tell()

	// 检测静音帧
	if tell >= totalBits {
		silence = 1
	} else if tell == 1 {
		silence = dec.DecBitLogp(15)
	} else {
		silence = 0
	}

	if silence != 0 {
		tell = totalBits
		dec.nbitsTotal += tell - dec.Tell()
	}

	// 解码后滤波器参数
	postfilterGain = 0
	postfilterPitch = 0
	postfilterTapset = 0
	if start == 0 && tell+16 <= totalBits {
		if dec.DecBitLogp(1) != 0 {
			octave := dec.DecUint(6)
			postfilterPitch = (16 << octave) + dec.DecBits(4+octave) - 1
			qg := dec.DecBits(3)
			if dec.Tell()+2 <= totalBits {
				postfilterTapset = dec.DecIcdf(tapsetIcdf)
			}
			postfilterGain = int(0.09375*float32(1<<15)) * (qg + 1)
		}
		tell = dec.Tell()
	}

	// 检测瞬态
	if LM > 0 && tell+3 <= totalBits {
		isTransient = dec.DecBitLogp(3)
		tell = dec.Tell()
	} else {
		isTransient = 0
	}

	if isTransient != 0 {
		shortBlocks = M
	} else {
		shortBlocks = 0
	}

	// 解码全局标志
	if tell+3 <= totalBits {
		intraEnerg = dec.DecBitLogp(3)
	} else {
		intraEnerg = 0
	}

	// 获取频带能量
	UnquantCoarseEnergy(mode, start, end, oldBandE, intraEnerg, dec, C, LM)

	// 解码时间-频率分辨率
	tfRes := make([]int, nbEBands)
	TfDecode(start, end, isTransient, tfRes, LM, dec)

	tell = dec.Tell()
	spreadDecision = SPREAD_NORMAL
	if tell+4 <= totalBits {
		spreadDecision = dec.DecIcdf(spreadIcdf)
	}

	// 初始化容量
	cap := make([]int, nbEBands)
	InitCaps(mode, cap, LM, C)

	// 动态分配
	offsets := make([]int, nbEBands)
	dynallocLogp = 6
	totalBits <<= BITRES
	tell = dec.TellFrac()
	for i = start; i < end; i++ {
		width := C * (eBands[i+1] - eBands[i]) << LM
		quanta := min(width<<BITRES, max(6<<BITRES, width))
		dynallocLoopLogp := dynallocLogp
		boost := 0
		for tell+(dynallocLoopLogp<<BITRES) < totalBits && boost < cap[i] {
			flag := dec.DecBitLogp(dynallocLoopLogp)
			tell = dec.TellFrac()
			if flag == 0 {
				break
			}
			boost += quanta
			totalBits -= quanta
			dynallocLoopLogp = 1
		}
		offsets[i] = boost
		if boost > 0 {
			dynallocLogp = max(2, dynallocLogp-1)
		}
	}

	// 精细量化
	fineQuant := make([]int, nbEBands)
	if tell+(6<<BITRES) <= totalBits {
		allocTrim = dec.DecIcdf(trimIcdf)
	} else {
		allocTrim = 5
	}

	bits = (len(data)*8)<<BITRES - dec.TellFrac() - 1
	if isTransient != 0 && LM >= 2 && bits >= (LM+2)<<BITRES {
		antiCollapseRsv = 1 << BITRES
	} else {
		antiCollapseRsv = 0
	}
	bits -= antiCollapseRsv

	// 计算分配
	pulses := make([]int, nbEBands)
	finePriority := make([]int, nbEBands)
	intensity = 0
	dualStereo = 0
	balance = 0
	codedBands := comm.ComputeAllocation(mode, start, end, offsets, cap, allocTrim,
		&intensity, &dualStereo, bits, &balance, pulses, fineQuant, finePriority,
		C, LM, dec)

	// 解码精细能量
	UnquantFineEnergy(mode, start, end, oldBandE, fineQuant, dec, C)

	// 移动解码内存
	for c := 0; c < CC; c++ {
		copy(d.decodeMem[c][:DECODE_BUFFER_SIZE-N+overlap/2], d.decodeMem[c][N:])
	}

	// 解码固定码本
	collapseMasks := make([]int16, C*nbEBands)
	X := make([][]int, C)
	for c := range X {
		X[c] = make([]int, N)
	}

	// 量化所有频带
	rng := d.rng
	QuantAllBands(0, mode, start, end, X[0], ternary(C == 2, X[1], nil), collapseMasks,
		nil, pulses, shortBlocks != 0, spreadDecision, dualStereo != 0, intensity, tfRes,
		len(data)*(8<<BITRES)-antiCollapseRsv, balance, dec, LM, codedBands, &rng)
	d.rng = rng

	// 反折叠
	if antiCollapseRsv > 0 {
		antiCollapseOn = dec.DecBits(1)
	}

	// 最终化能量
	UnquantEnergyFinalise(mode, start, end, oldBandE, fineQuant, finePriority,
		len(data)*8-dec.Tell(), dec, C)

	// 应用反折叠
	if antiCollapseOn != 0 {
		AntiCollapse(mode, X, collapseMasks, LM, C, N, start, end, oldBandE,
			oldLogE, oldLogE2, pulses, d.rng)
	}

	// 静音处理
	if silence != 0 {
		for i := 0; i < C*nbEBands; i++ {
			oldBandE[i] = -int(28.0 * float32(1<<DB_SHIFT))
		}
	}

	// 合成
	CeltSynthesis(mode, X, outSyn, outSynPtrs, oldBandE, start, effEnd,
		C, CC, isTransient != 0, LM, d.downsample, silence != 0)

	// 应用后滤波器
	for c := 0; c < CC; c++ {
		d.postfilterPeriod = max(d.postfilterPeriod, COMBFILTER_MINPERIOD)
		d.postfilterPeriodOld = max(d.postfilterPeriodOld, COMBFILTER_MINPERIOD)
		CombFilter(outSyn[c], outSynPtrs[c], outSyn[c], outSynPtrs[c],
			d.postfilterPeriodOld, d.postfilterPeriod, mode.shortMdctSize,
			d.postfilterGainOld, d.postfilterGain, d.postfilterTapsetOld, d.postfilterTapset,
			mode.window, overlap)

		if LM != 0 {
			CombFilter(outSyn[c], outSynPtrs[c]+mode.shortMdctSize,
				outSyn[c], outSynPtrs[c]+mode.shortMdctSize,
				d.postfilterPeriod, postfilterPitch, N-mode.shortMdctSize,
				d.postfilterGain, postfilterGain, d.postfilterTapset, postfilterTapset,
				mode.window, overlap)
		}
	}

	// 更新后滤波器状态
	d.postfilterPeriodOld = d.postfilterPeriod
	d.postfilterGainOld = d.postfilterGain
	d.postfilterTapsetOld = d.postfilterTapset
	d.postfilterPeriod = postfilterPitch
	d.postfilterGain = postfilterGain
	d.postfilterTapset = postfilterTapset

	// 单声道处理
	if C == 1 {
		copy(oldBandE[nbEBands:], oldBandE[:nbEBands])
	}

	// 更新背景噪声
	if isTransient == 0 {
		var maxBackgroundIncrease int
		if d.lossCount < 10 {
			maxBackgroundIncrease = M * int(0.001*float32(1<<DB_SHIFT))
		} else {
			maxBackgroundIncrease = int(1.0 * float32(1<<DB_SHIFT))
		}

		for i := 0; i < 2*nbEBands; i++ {
			backgroundLogE[i] = min(backgroundLogE[i]+maxBackgroundIncrease, oldBandE[i])
		}
	} else {
		for i := 0; i < 2*nbEBands; i++ {
			oldLogE[i] = min(oldLogE[i], oldBandE[i])
		}
	}

	// 清理频带
	for c := 0; c < 2; c++ {
		for i := 0; i < start; i++ {
			oldBandE[c*nbEBands+i] = 0
			val := -int(28.0 * float32(1<<DB_SHIFT))
			oldLogE[c*nbEBands+i] = val
			oldLogE2[c*nbEBands+i] = val
		}
		for i := end; i < nbEBands; i++ {
			oldBandE[c*nbEBands+i] = 0
			val := -int(28.0 * float32(1<<DB_SHIFT))
			oldLogE[c*nbEBands+i] = val
			oldLogE2[c*nbEBands+i] = val
		}
	}

	// 去加重
	Deemphasis(outSyn, outSynPtrs, pcm, N, CC, d.downsample, mode.preemph, d.preemphMemD[:], accum)
	d.lossCount = 0

	// 错误检查
	if dec.Tell() > 8*len(data) {
		return 0, errors.New("OPUS_INTERNAL_ERROR")
	}
	if dec.Error() != 0 {
		d.error = 1
	}

	return frameSize / d.downsample, nil
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

func ternary(cond bool, a, b []int) []int {
	if cond {
		return a
	}
	return b
}

func (d *CeltDecoder) SetStartBand(value int) error {
	if value < 0 || value >= d.mode.nbEBands {
		return errors.New("start band out of range")
	}
	d.start = value
	return nil
}

func (d *CeltDecoder) SetEndBand(value int) error {
	if value < 1 || value > d.mode.nbEBands {
		return errors.New("end band out of range")
	}
	d.end = value
	return nil
}

func (d *CeltDecoder) SetChannels(value int) error {
	if value < 1 || value > 2 {
		return errors.New("channel count must be 1 or 2")
	}
	d.streamChannels = value
	return nil
}

func (d *CeltDecoder) GetAndClearError() int {
	err := d.error
	d.error = 0
	return err
}

func (d *CeltDecoder) GetLookahead() int {
	if d.downsample == 0 {
		return 0
	}
	return d.overlap / d.downsample
}

func (d *CeltDecoder) GetPitch() int {
	return d.postfilterPeriod
}

func (d *CeltDecoder) GetMode() *CeltMode {
	return d.mode
}

func (d *CeltDecoder) SetSignalling(value int) {
	d.signalling = value
}

func (d *CeltDecoder) GetFinalRange() int {
	return d.rng
}
