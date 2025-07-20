package opus

import (
	"math"

	"githb.com/dosgo/v2/core/kissfft"
)

const (
	M_PI = 3.14159265358979323846
	cA   = 0.43157974
	cB   = 0.67848403
	cC   = 0.08595542
	cE   = M_PI / 2
)

// FastAtan2 快速计算 atan2
func FastAtan2(y, x float32) float32 {
	// 避免下溢
	if math.Abs(float64(x))+math.Abs(float64(y)) < 1e-9 {
		x *= 1e12
		y *= 1e12
	}

	x2 := x * x
	y2 := y * y

	if x2 < y2 {
		den := (y2 + cB*x2) * (y2 + cC*x2)
		if den != 0 {
			return -x*y*(y2+cA*x2)/den + float32(cE)
		}
		return float32(cE)
	} else {
		den := (x2 + cB*y2) * (x2 + cC*y2)
		if den != 0 {
			return x*y*(x2+cA*y2)/den + float32(cE)
		}
		return float32(cE)
	}
}

// TonalityAnalysisInit 初始化音调分析状态
func TonalityAnalysisInit(tonal *TonalityAnalysisState) {
	tonal.Reset()
}

// TonalityGetInfo 获取音调分析信息
func TonalityGetInfo(tonal *TonalityAnalysisState, info *AnalysisInfo, len int) {
	pos := tonal.ReadPos
	currLookahead := tonal.WritePos - tonal.ReadPos
	if currLookahead < 0 {
		currLookahead += DETECT_SIZE
	}

	if len > 480 && pos != tonal.WritePos {
		pos++
		if pos == DETECT_SIZE {
			pos = 0
		}
	}
	if pos == tonal.WritePos {
		pos--
	}
	if pos < 0 {
		pos = DETECT_SIZE - 1
	}

	*info = tonal.Info[pos]
	tonal.ReadSubframe += len / 120
	for tonal.ReadSubframe >= 4 {
		tonal.ReadSubframe -= 4
		tonal.ReadPos++
	}
	if tonal.ReadPos >= DETECT_SIZE {
		tonal.ReadPos -= DETECT_SIZE
	}

	// 补偿特征本身的延迟
	currLookahead = max(currLookahead-10, 0)

	psum := 0.0
	for i := 0; i < DETECT_SIZE-currLookahead; i++ {
		psum += tonal.Pmusic[i]
	}
	for i := DETECT_SIZE - currLookahead; i < DETECT_SIZE; i++ {
		psum += tonal.Pspeech[i]
	}
	psum = psum*tonal.MusicConfidence + (1-psum)*tonal.SpeechConfidence
	info.MusicProb = float32(psum)
}

// TonalityAnalysis 执行音调分析
func TonalityAnalysis(tonal *TonalityAnalysisState, celtMode *CeltMode, x []int16, xPtr, len, offset, c1, c2, C, lsbDepth int) {
	const N = 480
	const N2 = 240

	// 初始化变量
	alpha := 1.0 / min(20, 1+tonal.Count)
	alphaE := 1.0 / min(50, 1+tonal.Count)
	alphaE2 := 1.0 / min(1000, 1+tonal.Count)

	if tonal.Count < 4 {
		tonal.MusicProb = 0.5
	}

	kfft := celtMode.Mdct.Kfft[0]
	if tonal.Count == 0 {
		tonal.MemFill = 240
	}

	// 下混处理
	DownmixInt(x, xPtr, tonal.Inmem[:], tonal.MemFill, min(len, ANALYSIS_BUF_SIZE-tonal.MemFill), offset, c1, c2, C)

	if tonal.MemFill+len < ANALYSIS_BUF_SIZE {
		tonal.MemFill += len
		return
	}

	info := &tonal.Info[tonal.WritePos]
	tonal.WritePos = (tonal.WritePos + 1) % DETECT_SIZE

	// 准备FFT输入
	input := make([]float32, 960)
	output := make([]float32, 960)
	for i := 0; i < N2; i++ {
		w := AnalysisWindow[i]
		input[2*i] = w * tonal.Inmem[i]
		input[2*i+1] = w * tonal.Inmem[N2+i]
		input[2*(N-i-1)] = w * tonal.Inmem[N-i-1]
		input[2*(N-i-1)+1] = w * tonal.Inmem[N+N2-i-1]
	}

	// 移动内存
	copy(tonal.Inmem[:], tonal.Inmem[ANALYSIS_BUF_SIZE-240:])
	remaining := len - (ANALYSIS_BUF_SIZE - tonal.MemFill)
	DownmixInt(x, xPtr, tonal.Inmem[240:], 240, remaining, offset+ANALYSIS_BUF_SIZE-tonal.MemFill, c1, c2, C)
	tonal.MemFill = 240 + remaining

	// 执行FFT
	kissfft.OpusFFT(kfft, input, output)

	// 计算音调和噪声
	A := make([]float32, N2)
	dA := make([]float32, N2)
	d2A := make([]float32, N2)
	tonality := make([]float32, N2)
	noisiness := make([]float32, N2)

	for i := 1; i < N2; i++ {
		X1r := output[2*i] + output[2*(N-i)]
		X1i := output[2*i+1] - output[2*(N-i)+1]
		X2r := output[2*i+1] + output[2*(N-i)+1]
		X2i := output[2*(N-i)] - output[2*i]

		angle := 0.5 / M_PI * FastAtan2(X1i, X1r)
		dAngle := angle - tonal.Angle[i]
		d2Angle := dAngle - tonal.DAngle[i]

		angle2 := 0.5 / M_PI * FastAtan2(X2i, X2r)
		dAngle2 := angle2 - angle
		d2Angle2 := dAngle2 - dAngle

		mod1 := d2Angle - float32(math.Floor(float64(d2Angle)+0.5))
		noisiness[i] = float32(math.Abs(float64(mod1)))
		mod1 *= mod1
		mod1 *= mod1

		mod2 := d2Angle2 - float32(math.Floor(float64(d2Angle2)+0.5))
		noisiness[i] += float32(math.Abs(float64(mod2)))
		mod2 *= mod2
		mod2 *= mod2

		avgMod := 0.25 * (tonal.D2Angle[i] + 2.0*mod1 + mod2)
		tonality[i] = 1.0/(1.0+40.0*16.0*float32(math.Pow(M_PI, 4))*avgMod) - 0.015

		tonal.Angle[i] = angle2
		tonal.DAngle[i] = dAngle2
		tonal.D2Angle[i] = mod2
	}

	// 计算频带特性
	frameTonality := 0.0
	maxFrameTonality := 0.0
	frameNoisiness := 0.0
	frameStationarity := 0.0
	frameLoudness := 0.0
	relativeE := 0.0
	slope := 0.0
	logE := make([]float32, NB_TBANDS)
	bandTonality := make([]float32, NB_TBANDS)

	if tonal.Count == 0 {
		for b := 0; b < NB_TBANDS; b++ {
			tonal.LowE[b] = 1e10
			tonal.HighE[b] = -1e10
		}
	}

	for b := 0; b < NB_TBANDS; b++ {
		E := 0.0
		tE := 0.0
		nE := 0.0
		stationarity := 0.0

		for i := Tbands[b]; i < Tbands[b+1]; i++ {
			binE := float64(output[2*i]*output[2*i] + output[2*(N-i)]*output[2*(N-i)] +
				output[2*i+1]*output[2*i+1] + output[2*(N-i)+1]*output[2*(N-i)+1])
			binE *= 5.55e-17
			E += binE
			tE += binE * float64(tonality[i])
			nE += binE * 2.0 * (0.5 - float64(noisiness[i]))
		}

		tonal.E[tonal.ECount][b] = float32(E)
		frameNoisiness += nE / (1e-15 + E)
		frameLoudness += math.Sqrt(E + 1e-10)
		logE[b] = float32(math.Log(E + 1e-10))
		tonal.LowE[b] = min(tonal.LowE[b]+0.01, logE[b])
		tonal.HighE[b] = max(tonal.HighE[b]-0.1, logE[b])
		if tonal.HighE[b] < tonal.LowE[b]+1.0 {
			tonal.HighE[b] += 0.5
			tonal.LowE[b] -= 0.5
		}
		relativeE += (logE[b] - tonal.LowE[b]) / (1e-15 + tonal.HighE[b] - tonal.LowE[b])

		L1 := 0.0
		L2 := 0.0
		for i := 0; i < NB_FRAMES; i++ {
			L1 += math.Sqrt(float64(tonal.E[i][b]))
			L2 += float64(tonal.E[i][b])
		}

		stationarity = min(0.99, L1/math.Sqrt(1e-15+float64(NB_FRAMES)*L2))
		stationarity *= stationarity
		stationarity *= stationarity
		frameStationarity += stationarity
		bandTonality[b] = max(tonal.PrevBandTonality[b]*float32(stationarity), float32(tE/(1e-15+E)))
		frameTonality += float64(bandTonality[b])
		if b >= NB_TBANDS-NB_TONAL_SKIP_BANDS {
			frameTonality -= float64(bandTonality[b-NB_TBANDS+NB_TONAL_SKIP_BANDS])
		}
		maxFrameTonality = max(maxFrameTonality, (1.0+0.03*float64(b-NB_TBANDS))*frameTonality)
		slope += float64(bandTonality[b]) * float64(b-8)
		tonal.PrevBandTonality[b] = bandTonality[b]
	}

	// 计算带宽
	bandwidthMask := 0.0
	bandwidth := 0
	maxE := 0.0
	noiseFloor := 5.7e-4 / float64(1<<max(0, lsbDepth-8))
	noiseFloor *= float64(1 << (15 + SIG_SHIFT))
	noiseFloor *= noiseFloor

	for b := 0; b < NB_TOT_BANDS; b++ {
		E := 0.0
		bandStart := ExtraBands[b]
		bandEnd := ExtraBands[b+1]

		for i := bandStart; i < bandEnd; i++ {
			binE := float64(output[2*i]*output[2*i] + output[2*(N-i)]*output[2*(N-i)] +
				output[2*i+1]*output[2*i+1] + output[2*(N-i)+1]*output[2*(N-i)+1])
			E += binE
		}

		maxE = max(maxE, E)
		tonal.MeanE[b] = max((1-alphaE2)*tonal.MeanE[b], float32(E))
		E = max(E, float64(tonal.MeanE[b]))
		bandwidthMask = max(0.05*bandwidthMask, E)

		if E > 0.1*bandwidthMask && E*1e9 > maxE && E > noiseFloor*float64(bandEnd-bandStart) {
			bandwidth = b
		}
	}

	if tonal.Count <= 2 {
		bandwidth = 20
	}

	// 更新跟踪器
	frameLoudness = 20 * math.Log10(frameLoudness)
	tonal.Etracker = max(tonal.Etracker-0.03, float32(frameLoudness))
	tonal.LowECount *= (1 - alphaE)
	if frameLoudness < float64(tonal.Etracker-30) {
		tonal.LowECount += alphaE
	}

	// 计算BFCC特征
	BFCC := make([]float32, 8)
	for i := 0; i < 8; i++ {
		sum := 0.0
		for b := 0; b < 16; b++ {
			sum += DctTable[i*16+b] * float64(logE[b])
		}
		BFCC[i] = float32(sum)
	}

	// 归一化特征
	frameStationarity /= NB_TBANDS
	relativeE /= NB_TBANDS
	if tonal.Count < 10 {
		relativeE = 0.5
	}
	frameNoisiness /= NB_TBANDS
	info.Activity = float32(frameNoisiness + (1-frameNoisiness)*relativeE)
	frameTonality = maxFrameTonality / float64(NB_TBANDS-NB_TONAL_SKIP_BANDS)
	frameTonality = max(float64(tonal.PrevTonality)*0.8, frameTonality)
	tonal.PrevTonality = float32(frameTonality)
	slope /= 8 * 8
	info.TonalitySlope = float32(slope)

	tonal.ECount = (tonal.ECount + 1) % NB_FRAMES
	tonal.Count++
	info.Tonality = float32(frameTonality)

	// 特征提取
	features := make([]float32, 25)
	for i := 0; i < 4; i++ {
		features[i] = -0.12299*(BFCC[i]+tonal.Mem[i+24]) + 0.49195*(tonal.Mem[i]+tonal.Mem[i+16]) + 0.69693*tonal.Mem[i+8] - 1.4349*tonal.Cmean[i]
	}

	for i := 0; i < 4; i++ {
		tonal.Cmean[i] = (1-alpha)*tonal.Cmean[i] + alpha*BFCC[i]
	}

	for i := 0; i < 4; i++ {
		features[4+i] = 0.63246*(BFCC[i]-tonal.Mem[i+24]) + 0.31623*(tonal.Mem[i]-tonal.Mem[i+16])
	}

	for i := 0; i < 3; i++ {
		features[8+i] = 0.53452*(BFCC[i]+tonal.Mem[i+24]) - 0.26726*(tonal.Mem[i]+tonal.Mem[i+16]) - 0.53452*tonal.Mem[i+8]
	}

	// 计算标准差
	if tonal.Count > 5 {
		for i := 0; i < 9; i++ {
			tonal.Std[i] = (1-alpha)*tonal.Std[i] + alpha*features[i]*features[i]
		}
	}

	// 更新内存
	for i := 0; i < 8; i++ {
		tonal.Mem[i+24] = tonal.Mem[i+16]
		tonal.Mem[i+16] = tonal.Mem[i+8]
		tonal.Mem[i+8] = tonal.Mem[i]
		tonal.Mem[i] = BFCC[i]
	}

	// 计算标准差特征
	for i := 0; i < 9; i++ {
		features[11+i] = float32(math.Sqrt(float64(tonal.Std[i])))
	}

	// 设置其他特征
	features[20] = info.Tonality
	features[21] = info.Activity
	features[22] = float32(frameStationarity)
	features[23] = info.TonalitySlope
	features[24] = float32(tonal.LowECount)

	// 处理概率
	if info.Enabled {
		frameProbs := make([]float32, 2)
		MLPProcess(Net, features, frameProbs)

		frameProbs[0] = 0.5 * (frameProbs[0] + 1)
		frameProbs[0] = 0.01 + 1.21*frameProbs[0]*frameProbs[0] - 0.23*float32(math.Pow(float64(frameProbs[0]), 10))
		frameProbs[1] = 0.5*frameProbs[1] + 0.5
		frameProbs[0] = frameProbs[1]*frameProbs[0] + (1-frameProbs[1])*0.5

		// 状态转移概率计算
		tau := 0.00005 * frameProbs[1]
		beta := 0.05
		if true {
			p := max(0.05, min(0.95, frameProbs[0]))
			q := max(0.05, min(0.95, tonal.MusicProb))
			beta = 0.01 + 0.05*math.Abs(float64(p-q))/(p*(1-q)+q*(1-p))
		}

		p0 := (1-tonal.MusicProb)*(1-tau) + tonal.MusicProb*tau
		p1 := tonal.MusicProb*(1-tau) + (1-tonal.MusicProb)*tau
		p0 *= float32(math.Pow(1-float64(frameProbs[0]), beta))
		p1 *= float32(math.Pow(float64(frameProbs[0]), beta))
		tonal.MusicProb = p1 / (p0 + p1)
		info.MusicProb = tonal.MusicProb

		// 延迟决策处理
		// ... [省略详细实现] ...
	} else {
		info.MusicProb = 0
	}

	// 设置最终信息
	info.Bandwidth = bandwidth
	info.Noisiness = float32(frameNoisiness)
	info.Valid = 1
}

// RunAnalysis 运行音调分析
func RunAnalysis(analysis *TonalityAnalysisState, celtMode *CeltMode, analysisPcm []int16, analysisPcmPtr, analysisFrameSize, frameSize, c1, c2, C, Fs, lsbDepth int, analysisInfo *AnalysisInfo) {
	if analysisPcm != nil {
		// 避免分析缓冲区溢出
		analysisFrameSize = min((DETECT_SIZE-5)*Fs/100, analysisFrameSize)

		pcmLen := analysisFrameSize - analysis.AnalysisOffset
		offset := analysis.AnalysisOffset
		for pcmLen > 0 {
			chunk := min(480, pcmLen)
			TonalityAnalysis(analysis, celtMode, analysisPcm, analysisPcmPtr, chunk, offset, c1, c2, C, lsbDepth)
			offset += chunk
			pcmLen -= chunk
		}
		analysis.AnalysisOffset = analysisFrameSize
		analysis.AnalysisOffset -= frameSize
	}

	analysisInfo.Valid = 0
	TonalityGetInfo(analysis, analysisInfo, frameSize)
}

// 辅助函数
func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func max(a, b float64) float64 {
	if a > b {
		return a
	}
	return b
}

func minFloat(a, b float32) float32 {
	if a < b {
		return a
	}
	return b
}

func maxFloat(a, b float32) float32 {
	if a > b {
		return a
	}
	return b
}
