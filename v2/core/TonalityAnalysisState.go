package opus

// TonalityAnalysisState 音调分析状态结构体
type TonalityAnalysisState struct {
	Enabled               bool
	Angle                 [240]float32
	DAngle                [240]float32
	D2Angle               [240]float32
	Inmem                 [ANALYSIS_BUF_SIZE]int32
	MemFill               int32
	PrevBandTonality      [NB_TBANDS]float32
	PrevTonality          float32
	E                     [NB_FRAMES][NB_TBANDS]float32
	LowE                  [NB_TBANDS]float32
	HighE                 [NB_TBANDS]float32
	MeanE                 [NB_TOT_BANDS]float32
	Mem                   [32]float32
	Cmean                 [8]float32
	Std                   [9]float32
	MusicProb             float32
	Etracker              float32
	LowECount             float32
	ECount                int32
	LastMusic             int32
	LastTransition        int32
	Count                 int32
	SubframeMem           [3]float32
	AnalysisOffset        int32
	Pspeech               [DETECT_SIZE]float32
	Pmusic                [DETECT_SIZE]float32
	SpeechConfidence      float32
	MusicConfidence       float32
	SpeechConfidenceCount int32
	MusicConfidenceCount  int32
	WritePos              int32
	ReadPos               int32
	ReadSubframe          int32
	Info                  [DETECT_SIZE]AnalysisInfo
}

// NewTonalityAnalysisState 创建并初始化状态实例
func NewTonalityAnalysisState() *TonalityAnalysisState {
	state := &TonalityAnalysisState{
		// Go 数组自动初始化为零值，无需显式 memset
	}
	for i := range state.Info {
		state.Info[i] = AnalysisInfo{}
	}
	return state
}

// Reset 重置所有状态到初始值
func (s *TonalityAnalysisState) Reset() {
	s.Enabled = false
	s.Angle = [240]float32{}
	s.DAngle = [240]float32{}
	s.D2Angle = [240]float32{}
	s.Inmem = [ANALYSIS_BUF_SIZE]int32{}
	s.MemFill = 0
	s.PrevBandTonality = [NB_TBANDS]float32{}
	s.PrevTonality = 0
	for i := range s.E {
		s.E[i] = [NB_TBANDS]float32{}
	}
	s.LowE = [NB_TBANDS]float32{}
	s.HighE = [NB_TBANDS]float32{}
	s.MeanE = [NB_TOT_BANDS]float32{}
	s.Mem = [32]float32{}
	s.Cmean = [8]float32{}
	s.Std = [9]float32{}
	s.MusicProb = 0
	s.Etracker = 0
	s.LowECount = 0
	s.ECount = 0
	s.LastMusic = 0
	s.LastTransition = 0
	s.Count = 0
	s.SubframeMem = [3]float32{}
	s.AnalysisOffset = 0
	s.Pspeech = [DETECT_SIZE]float32{}
	s.Pmusic = [DETECT_SIZE]float32{}
	s.SpeechConfidence = 0
	s.MusicConfidence = 0
	s.SpeechConfidenceCount = 0
	s.MusicConfidenceCount = 0
	s.WritePos = 0
	s.ReadPos = 0
	s.ReadSubframe = 0
	for i := range s.Info {
		s.Info[i].Reset()
	}
}
