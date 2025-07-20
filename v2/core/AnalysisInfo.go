package opus

// AnalysisInfo 存储音频分析结果的结构体
type AnalysisInfo struct {
	Enabled       bool    // 是否启用分析
	Valid         int     // 分析结果的有效性标志
	Tonality      float32 // 音调性强度 (0-1)
	TonalitySlope float32 // 音调性变化斜率
	Noisiness     float32 // 噪声强度 (0-1)
	Activity      float32 // 音频活动强度 (0-1)
	MusicProb     float32 // 音乐信号概率 (0-1)
	Bandwidth     int     // 检测到的音频带宽
}

// NewAnalysisInfo 创建并初始化AnalysisInfo实例
func NewAnalysisInfo() *AnalysisInfo {
	return &AnalysisInfo{
		Enabled:   false,
		Bandwidth: -1, // 初始化为无效值
	}
}

// Assign 从另一个实例复制分析数据
func (a *AnalysisInfo) Assign(other *AnalysisInfo) {
	if other == nil {
		return
	}

	a.Valid = other.Valid
	a.Tonality = other.Tonality
	a.TonalitySlope = other.TonalitySlope
	a.Noisiness = other.Noisiness
	a.Activity = other.Activity
	a.MusicProb = other.MusicProb
	a.Bandwidth = other.Bandwidth
	// 注意: enabled 状态不被复制
}

// Reset 重置分析数据到初始状态
func (a *AnalysisInfo) Reset() {
	a.Valid = 0
	a.Tonality = 0
	a.TonalitySlope = 0
	a.Noisiness = 0
	a.Activity = 0
	a.MusicProb = 0
	a.Bandwidth = -1 // 重置为无效值
}
