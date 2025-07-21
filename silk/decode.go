package silk

import (
	"errors"
	"math"
	"math/bits"
	"sort"

	"github.com/dosgo/goOpus/comm"
)

// 带宽类型
type Bandwidth int

const (
	BandwidthNarrow Bandwidth = iota
	BandwidthMedium
	BandwidthWide
	BandwidthSuperwide
	BandwidthFull
)

// SilkInfo 保存配置信息
type SilkInfo struct {
	Bandwidth Bandwidth
	Subframes int
	SfSize    int
	FSize     int

	Weight0 float32
	Weight1 float32
	Prev0   float32
	Prev1   float32
}

// Silk 处理主要解码过程
type Silk struct {
	Stereo      bool
	StereoOut   bool
	Frames      int
	FrameLen    int
	SubframeLen int
	Info        SilkInfo

	MidFrame    *SilkFrame
	SideFrame   *SilkFrame
	LeftOutbuf  []float32
	RightOutbuf []float32
}

// SubFrame 包含音频子帧数据
type SubFrame struct {
	Gain     float32
	PitchLag int32
	LtpTaps  [5]float32
}

// FrameType 表示编码特性
type FrameType struct {
	Active bool
	Voiced bool
	High   bool
}

func (ft FrameType) VoicedIndex() int {
	if ft.Voiced {
		return 1
	}
	return 0
}

func (ft FrameType) SignalTypeIndex() int {
	idx := 0
	if ft.Voiced {
		idx++
	}
	if ft.Active {
		idx++
	}
	return idx
}

func (ft FrameType) QoffsetTypeIndex() int {
	if ft.High {
		return 1
	}
	return 0
}

// Band 表示频段信息
type Band struct {
	Order           int
	Step            int32
	Stage1          []*comm.ICDFContext
	Map             [][]*comm.ICDFContext
	PredWeight      [][]uint8
	PredWeightIndex [][]int
	Weight          [][]uint16
	Codebook        [][]uint8
	MinSpacing      []int16
	Ordering        []uint8
}

func log2Lin(val int) int {
	i := 1 << (val >> 7)
	f := val & 127
	return i + ((-174*f*(128-f)>>16)+f)*(i>>7)
}

func mulShift32(a, b int32, bits uint) int32 {
	return int32((int64(a) * int64(b)) >> bits)
}

func mulRound32(a, b int32, bits uint) int32 {
	product := int64(a) * int64(b)
	half := int64(1) << (bits - 1)
	return int32((product + half) >> bits)
}

func (b *Band) Stabilize(nlsfs []int16) {
	n := b.Order
	for iter := 0; iter < 20; iter++ {
		minDiff := int32(0)
		minK := -1

		for i := 0; i <= n; i++ {
			var low, high int32
			if i > 0 {
				low = int32(nlsfs[i-1])
			}
			if i < n {
				high = int32(nlsfs[i])
			} else {
				high = 32768
			}

			diff := high - low - int32(b.MinSpacing[i])
			if minK < 0 || diff < minDiff {
				minDiff = diff
				minK = i
			}
		}

		if minDiff == 0 {
			return
		}

		switch {
		case minK == 0:
			nlsfs[0] = b.MinSpacing[0]
		case minK == n:
			nlsfs[n-1] = int16(32768 - int32(b.MinSpacing[n]))
		default:
			halfDelta := int32(b.MinSpacing[minK]) / 2

			minCenter := int32(0)
			for j := 0; j < minK; j++ {
				minCenter += int32(b.MinSpacing[j])
			}
			minCenter += halfDelta

			maxCenter := int32(32768)
			for j := minK + 1; j <= n; j++ {
				maxCenter -= int32(b.MinSpacing[j])
			}
			maxCenter -= halfDelta

			// 计算中心点
			sum := int32(nlsfs[minK-1]) + int32(nlsfs[minK])
			center := sum/2 + sum%2

			// 钳位中心点
			if center < minCenter {
				center = minCenter
			} else if center > maxCenter {
				center = maxCenter
			}

			nlsfs[minK-1] = int16(center - halfDelta)
			nlsfs[minK] = int16(center + halfDelta)
		}
	}

	// 迭代后的最终稳定处理
	sort.Slice(nlsfs[:n], func(i, j int) bool {
		return nlsfs[i] < nlsfs[j]
	})

	prev := int32(0)
	for i := 0; i < n; i++ {
		minVal := int32(prev + int32(b.MinSpacing[i]))
		if int32(nlsfs[i]) < minVal {
			nlsfs[i] = int16(minVal)
		}
		prev = int32(nlsfs[i])
	}

	next := int32(32768)
	for i := n - 1; i >= 0; i-- {
		minVal := next - int32(b.MinSpacing[i+1])
		if int32(nlsfs[i]) > minVal {
			nlsfs[i] = int16(minVal)
		}
		next = int32(nlsfs[i])
	}
}

// is_stable 检查LPC系数是否稳定
func (b *Band) IsStable(lpcs []int16) bool {
	dc_resp := 0
	even := make([]int32, b.Order)
	odd := make([]int32, b.Order)
	invgain := int(1 << 30)

	// 计算直流响应
	for i, lpc := range lpcs {
		l := int32(lpc)
		dc_resp += int(l)
		even[i] = l * 4096
	}

	if dc_resp > 4096 {
		return false
	}

	k := b.Order - 1
	a := even[k]

	for {
		if abs32(a) > 16773022 {
			return false
		}

		rc := -a * 128
		div := (1 << 30) - mulShift32(rc, rc, 32)

		invgain = int(mulShift32(int32(invgain>>32), div, 32) << 34)

		if k == 0 {
			return invgain >= 107374
		}

		b1 := celtIlog2(div)
		b2 := b1 - 16
		divShifted := div >> (b2 + 1)
		inv := int32(((1 << 29) - 1) / divShifted)
		err := (1 << 29) - mulShift32(div<<(15-b2), inv, 16)
		gain := (inv << 16) + mulShift32(err, inv, 13)

		prev := even
		cur := odd
		if k&1 != 0 {
			prev = odd
			cur = even
		}

		for j := 0; j < k; j++ {
			v := prev[j] - mulShift32(prev[k-j-1], rc, 31)
			cur[j] = mulShift32(v, gain, uint(b1))
		}

		k--
		a = cur[k]
	}
}

// range_limit 限制LPC系数的范围
func (b *Band) RangeLimit(lpcs []float32, a []int32) {
	lpc := make([]int16, b.Order)
	deadline := true

	for iter := 0; iter < 10; iter++ {
		maxabs := int32(0)
		k := 0

		// 找到最大绝对值及其索引
		for i, val := range a {
			absVal := abs32(val)
			if absVal >= maxabs {
				maxabs = absVal
				k = i
			}
		}

		maxabs = (maxabs + 16) >> 5
		if maxabs > 32767 {
			maxVal := int32(math.Max(float64(maxabs), 163838))
			start := 65470 - ((maxVal-32767)<<14)/(maxVal*int32(k+1)>>2)
			chirp := start

			for i := range a {
				a[i] = mulRound32(a[i], chirp, 16)
				chirp = (start*chirp + 32768) >> 16
			}
		} else {
			deadline = false
			break
		}
	}

	if deadline {
		for i, val := range a {
			v16 := (val + 16) >> 5
			if v16 > math.MaxInt16 {
				v16 = math.MaxInt16
			} else if v16 < math.MinInt16 {
				v16 = math.MinInt16
			}
			lpc[i] = int16(v16)
			a[i] = v16 << 5
		}
	} else {
		for i, val := range a {
			lpc[i] = int16((val + 16) >> 5)
		}
	}

	for i := 1; i <= 16; i++ {
		if b.IsStable(lpc) {
			break
		}

		start := uint32(65536 - (1 << i))
		chirp := start

		for j := range a {
			a[j] = mulRound32(a[j], int32(chirp), 16)
			lpc[j] = int16((a[j] + 16) >> 5)
			chirp = (start*chirp + 32768) >> 16
		}
	}

	for i := range lpcs {
		lpcs[i] = float32(lpc[i]) / 4096.0
	}
}

// lsf_to_lpc 将线谱频率转换为线性预测系数
func (b *Band) LsfToLpc(lpcs []float32, nlsfs []int16) {
	lsps := make([]int32, b.Order)
	p := make([]int32, b.Order/2+1)
	q := make([]int32, b.Order/2+1)

	for i, ord := range b.Ordering {
		nlsf := nlsfs[i]
		idx := nlsf >> 8
		off := nlsf & 255

		cos := int32(COSINE[idx])
		nextCos := int32(COSINE[idx+1])

		lsps[ord] = (cos*256 + (nextCos-cos)*int32(off) + 4) >> 3
	}

	p[0] = 65536
	q[0] = 65536
	p[1] = -lsps[0]
	q[1] = -lsps[1]

	for i := 0; i < b.Order/2-1; i++ {
		lsp0 := lsps[2+i*2]
		lsp1 := lsps[3+i*2]

		p[i+2] = p[i]*2 - mulRound32(lsp0, p[i+1], 16)
		q[i+2] = q[i]*2 - mulRound32(lsp1, q[i+1], 16)

		for j := i + 1; j > 0; j-- {
			p[j] += p[j-2] - mulRound32(lsp0, p[j-1], 16)
			q[j] += q[j-2] - mulRound32(lsp1, q[j-1], 16)
		}

		p[1] -= lsp0
		q[1] -= lsp1
	}

	a := make([]int32, b.Order)
	halfOrder := b.Order / 2
	for i := 0; i < halfOrder; i++ {
		ps := p[i] + p[i+1]
		qs := q[i+1] - q[i]
		a[i] = -ps - qs
		a[b.Order-i-1] = -ps + qs
	}

	b.RangeLimit(lpcs, a)
}

// 辅助函数
func abs32(x int32) int32 {
	if x < 0 {
		return -x
	}
	return x
}

// celtIlog2 计算整数的以2为底的对数
func celtIlog2(x int32) int {
	if x <= 0 {
		return 0
	}
	return 31 - bits.LeadingZeros32(uint32(x))
}

// PitchLag 接口定义了音高延迟相关功能
type PitchLag interface {
	LowPart() *comm.ICDFContext
	MinLag() uint16
	MaxLag() uint16
	Scale() uint16
	Offset() [][][]int8
	Contour() []*comm.ICDFContext
}
type CommPitchLag struct {
	LOW_PART *comm.ICDFContext
	MIN_LAG  uint16
	MAX_LAG  uint16
	SCALE    uint16
	OFFSET   [][][]int8
	CONTOUR  []*comm.ICDFContext
}

func (pitchLag CommPitchLag) LowPart() *comm.ICDFContext {
	return pitchLag.LOW_PART
}
func (pitchLag CommPitchLag) MinLag() uint16 {
	return pitchLag.MIN_LAG
}
func (pitchLag CommPitchLag) MaxLag() uint16 {
	return pitchLag.MAX_LAG
}
func (pitchLag CommPitchLag) Scale() uint16 {
	return pitchLag.SCALE
}
func (pitchLag CommPitchLag) Offset() [][][]int8 {
	return pitchLag.OFFSET
}
func (pitchLag CommPitchLag) Contour() []*comm.ICDFContext {
	return pitchLag.CONTOUR
}

var NB_MB = Band{
	Order:           10,
	Step:            11796,
	Stage1:          LSF_STAGE1_NB_MB,
	Map:             LSF_STAGE2_NB_MB.MAP,
	PredWeight:      LSF_PRED_WEIGHT_NB_MB,
	PredWeightIndex: LSF_PRED_WEIGHT_INDEX_NB_MB,
	Weight:          LSF_WEIGHT_NB_MB,
	Codebook:        LSF_CODEBOOK_NB_MB,
	MinSpacing:      LSF_MIN_SPACING_NB_MB,
	Ordering:        LSF_ORDERING_NB_MB,
}

var WB = Band{
	Order:           16,
	Step:            9830,
	Stage1:          LSF_STAGE1_WB,
	Map:             lsf_stage2_wb.MAP,
	PredWeight:      LSF_PRED_WEIGHT_WB,
	PredWeightIndex: LSF_PRED_WEIGHT_INDEX_WB,
	Weight:          LSF_WEIGHT_WB,
	Codebook:        LSF_CODEBOOK_WB,
	MinSpacing:      LSF_MIN_SPACING_WB,
	Ordering:        LSF_ORDERING_WB,
}

// 音高常量定义
var PITCH_HIGH_PART = &comm.ICDFContext{
	Total: 256,
	Dist: []int{
		3, 6, 12, 23, 44, 74, 106, 125, 136, 146, 158, 171, 184, 196, 207, 216, 224, 231, 237, 241,
		243, 245, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256,
	},
}

var PITCH_OFFSET_NB = [][][]int8{
	{
		{0, 0}, {1, 0}, {0, 1},
	},
	{
		{0, 0, 0, 0},
		{2, 1, 0, -1},
		{-1, 0, 1, 2},
		{-1, 0, 0, 1},
		{-1, 0, 0, 0},
		{0, 0, 0, 1},
		{0, 0, 1, 1},
		{1, 1, 0, 0},
		{1, 0, 0, 0},
		{0, 0, 0, -1},
		{1, 0, 0, -1},
	},
}

var PITCH_CONTOUR_NB = []*comm.ICDFContext{
	{
		Total: 256,
		Dist:  []int{143, 193, 256},
	},
	{
		Total: 256,
		Dist:  []int{68, 80, 101, 118, 137, 159, 189, 213, 230, 246, 256},
	},
}

var PITCH_OFFSET_MB_WB = [][][]int8{
	{
		{0, 0},
		{0, 1},
		{1, 0},
		{-1, 1},
		{1, -1},
		{-1, 2},
		{2, -1},
		{-2, 2},
		{2, -2},
		{-2, 3},
		{3, -2},
		{-3, 3},
	},
	{
		{0, 0, 0, 0},
		{0, 0, 1, 1},
		{1, 1, 0, 0},
		{-1, 0, 0, 0},
		{0, 0, 0, 1},
		{1, 0, 0, 0},
		{-1, 0, 0, 1},
		{0, 0, 0, -1},
		{-1, 0, 1, 2},
		{1, 0, 0, -1},
		{-2, -1, 1, 2},
		{2, 1, 0, -1},
		{-2, 0, 0, 2},
		{-2, 0, 1, 3},
		{2, 1, -1, -2},
		{-3, -1, 1, 3},
		{2, 0, 0, -2},
		{3, 1, 0, -2},
		{-3, -1, 2, 4},
		{-4, -1, 1, 4},
		{3, 1, -1, -3},
		{-4, -1, 2, 5},
		{4, 2, -1, -3},
		{4, 1, -1, -4},
		{-5, -1, 2, 6},
		{5, 2, -1, -4},
		{-6, -2, 2, 6},
		{-5, -2, 2, 5},
		{6, 2, -1, -5},
		{-7, -2, 3, 8},
		{6, 2, -2, -6},
		{5, 2, -2, -5},
		{8, 3, -2, -7},
		{-9, -3, 3, 9},
	},
}

var PITCH_CONTOUR_MB_WB = []*comm.ICDFContext{
	{
		Total: 256,
		Dist:  []int{91, 137, 176, 195, 209, 221, 229, 236, 242, 247, 252, 256},
	},
	{
		Total: 256,
		Dist: []int{
			33, 55, 73, 89, 104, 118, 132, 145, 158, 168, 177, 186, 194, 200, 206, 212, 217, 221,
			225, 229, 232, 235, 238, 240, 242, 244, 246, 248, 250, 252, 253, 254, 255, 256,
		},
	},
}

// NB 带宽实现 PitchLag 接口
var NB_PitchLag = &CommPitchLag{
	LOW_PART: &comm.ICDFContext{
		Total: 256,
		Dist:  []int{64, 128, 192, 256},
	},
	MIN_LAG: 16,
	MAX_LAG: 144,
	SCALE:   4,
	OFFSET:  PITCH_OFFSET_NB,
	CONTOUR: PITCH_CONTOUR_NB,
}
var MB_PitchLag = &CommPitchLag{
	LOW_PART: &comm.ICDFContext{
		Total: 256,
		Dist:  []int{43, 85, 128, 171, 213, 256},
	},
	MIN_LAG: 24,
	MAX_LAG: 216,
	SCALE:   6,
	OFFSET:  PITCH_OFFSET_MB_WB,
	CONTOUR: PITCH_CONTOUR_MB_WB,
}

var WB_PitchLag = &CommPitchLag{
	LOW_PART: &comm.ICDFContext{
		Total: 256,
		Dist:  []int{32, 64, 96, 128, 160, 192, 224, 256},
	},
	MIN_LAG: 32,
	MAX_LAG: 288,
	SCALE:   8,
	OFFSET:  PITCH_OFFSET_MB_WB,
	CONTOUR: PITCH_CONTOUR_MB_WB,
}

// LTP_PERIODICITY 是长期周期性参数
var LTP_PERIODICITY = &comm.ICDFContext{
	Total: 256,
	Dist:  []int{77, 157, 256},
}

// LTP_FILTER 长时预测滤波器参数
var LTP_FILTER = []*comm.ICDFContext{
	{
		Total: 256,
		Dist:  []int{185, 200, 213, 226, 235, 244, 250, 256},
	},
	{
		Total: 256,
		Dist:  []int{57, 91, 112, 132, 147, 160, 172, 185, 195, 205, 214, 224, 233, 241, 248, 256},
	},
	{
		Total: 256,
		Dist: []int{
			15, 31, 45, 57, 69, 81, 92, 103, 114, 124, 133, 142, 151, 160, 168, 176, 184, 192, 199,
			206, 212, 218, 223, 227, 232, 236, 240, 244, 247, 251, 254, 256,
		},
	},
}

// LTP_TAPS 长时预测系数
var LTP_TAPS = [][][]int8{
	{
		{4, 6, 24, 7, 5},
		{0, 0, 2, 0, 0},
		{12, 28, 41, 13, -4},
		{-9, 15, 42, 25, 14},
		{1, -2, 62, 41, -9},
		{-10, 37, 65, -4, 3},
		{-6, 4, 66, 7, -8},
		{16, 14, 38, -3, 33},
	},
	{
		{13, 22, 39, 23, 12},
		{-1, 36, 64, 27, -6},
		{-7, 10, 55, 43, 17},
		{1, 1, 8, 1, 1},
		{6, -11, 74, 53, -9},
		{-12, 55, 76, -12, 8},
		{-3, 3, 93, 27, -4},
		{26, 39, 59, 3, -8},
		{2, 0, 77, 11, 9},
		{-8, 22, 44, -6, 7},
		{40, 9, 26, 3, 9},
		{-7, 20, 101, -7, 4},
		{3, -8, 42, 26, 0},
		{-15, 33, 68, 2, 23},
		{-2, 55, 46, -2, 15},
		{3, -1, 21, 16, 41},
	},
	{
		{-6, 27, 61, 39, 5},
		{-11, 42, 88, 4, 1},
		{-2, 60, 65, 6, -4},
		{-1, -5, 73, 56, 1},
		{-9, 19, 94, 29, -9},
		{0, 12, 99, 6, 4},
		{8, -19, 102, 46, -13},
		{3, 2, 13, 3, 2},
		{9, -21, 84, 72, -18},
		{-11, 46, 104, -22, 8},
		{18, 38, 48, 23, 0},
		{-16, 70, 83, -21, 11},
		{5, -11, 117, 22, -8},
		{-6, 23, 117, -12, 3},
		{3, -8, 95, 28, 4},
		{-10, 15, 77, 60, -15},
		{-1, 4, 124, 2, -4},
		{3, 38, 84, 24, -25},
		{2, 13, 42, 13, 31},
		{21, -4, 56, 46, -1},
		{-1, 35, 79, -13, 19},
		{-7, 65, 88, -9, -14},
		{20, 4, 81, 49, -29},
		{20, 0, 75, 3, -17},
		{5, -9, 44, 92, -8},
		{1, -3, 22, 69, 31},
		{-6, 95, 41, -12, 5},
		{39, 67, 16, -4, 1},
		{0, -6, 120, 55, -36},
		{-13, 44, 122, 4, -24},
		{81, 5, 11, 3, 7},
		{2, 0, 9, 10, 88},
	},
}

// LTP_SCALE 长时预测缩放系数
var LTP_SCALE = []uint16{15565, 12288, 8192}

// LTP_SCALE_INDEX 长时预测缩放索引
var LTP_SCALE_INDEX = &comm.ICDFContext{
	Total: 256,
	Dist:  []int{128, 192, 256},
}

// LTP_ORDER 长时预测阶数
const LTP_ORDER = 5

// RES_HISTORY 残差历史长度
const RES_HISTORY = 288 + LTP_ORDER/2

// LPC_HISTORY 线性预测历史长度
const LPC_HISTORY = 322

// LCG_SEED 线性同余生成器种子
var LCG_SEED = &comm.ICDFContext{
	Total: 256,
	Dist:  []int{64, 128, 192, 256},
}

// ShellBlock 接口定义了码本块相关功能
type ShellBlock struct {
	ShellBlocks []uint8
}

var NBShellBlock = ShellBlock{[]uint8{5, 10}}

var MBShellBlock = ShellBlock{[]uint8{8, 15}}

var WBShellBlock = ShellBlock{[]uint8{10, 20}}

// EXC_RATE 激励率参数
var EXC_RATE = []*comm.ICDFContext{
	{
		Total: 256,
		Dist:  []int{15, 66, 78, 124, 169, 182, 215, 242, 256},
	},
	{
		Total: 256,
		Dist:  []int{33, 63, 99, 116, 150, 199, 217, 238, 256},
	},
}

// PULSE_COUNT 脉冲计数参数
var PULSE_COUNT = []*comm.ICDFContext{
	{
		Total: 256,
		Dist: []int{
			131, 205, 230, 238, 241, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256,
		},
	},
	{
		Total: 256,
		Dist: []int{
			58, 151, 211, 234, 241, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256,
		},
	},
	{
		Total: 256,
		Dist: []int{
			43, 94, 140, 173, 197, 213, 224, 232, 238, 241, 244, 247, 249, 250, 251, 253, 254, 256,
		},
	},
	{
		Total: 256,
		Dist: []int{
			17, 69, 140, 197, 228, 240, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256,
		},
	},
	{
		Total: 256,
		Dist: []int{
			6, 27, 68, 121, 170, 205, 226, 237, 243, 246, 248, 250, 251, 252, 253, 254, 255, 256,
		},
	},
	{
		Total: 256,
		Dist: []int{
			7, 21, 43, 71, 100, 128, 153, 173, 190, 203, 214, 223, 230, 235, 239, 243, 246, 256,
		},
	},
	{
		Total: 256,
		Dist: []int{
			2, 7, 21, 50, 92, 138, 179, 210, 229, 240, 246, 249, 251, 252, 253, 254, 255, 256,
		},
	},
	{
		Total: 256,
		Dist: []int{
			1, 3, 7, 17, 36, 65, 100, 137, 171, 199, 219, 233, 241, 246, 250, 252, 254, 256,
		},
	},
	{
		Total: 256,
		Dist: []int{
			1, 3, 5, 10, 19, 33, 53, 77, 104, 132, 158, 181, 201, 216, 227, 235, 241, 256,
		},
	},
	{
		Total: 256,
		Dist: []int{
			1, 2, 3, 9, 36, 94, 150, 189, 214, 228, 238, 244, 247, 250, 252, 253, 254, 256,
		},
	},
	{
		Total: 256,
		Dist: []int{
			2, 3, 9, 36, 94, 150, 189, 214, 228, 238, 244, 247, 250, 252, 253, 254, 256, 256,
		},
	},
}

// PULSE_LOCATION 脉冲位置分布参数
var PULSE_LOCATION = [][]*comm.ICDFContext{
	{
		{
			Total: 256,
			Dist:  []int{126, 256},
		},
		{
			Total: 256,
			Dist:  []int{56, 198, 256},
		},
		{
			Total: 256,
			Dist:  []int{25, 126, 230, 256},
		},
		{
			Total: 256,
			Dist:  []int{12, 72, 180, 244, 256},
		},
		{
			Total: 256,
			Dist:  []int{7, 42, 126, 213, 250, 256},
		},
		{
			Total: 256,
			Dist:  []int{4, 24, 83, 169, 232, 253, 256},
		},
		{
			Total: 256,
			Dist:  []int{3, 15, 53, 125, 200, 242, 254, 256},
		},
		{
			Total: 256,
			Dist:  []int{2, 10, 35, 89, 162, 221, 248, 255, 256},
		},
		{
			Total: 256,
			Dist:  []int{2, 7, 24, 63, 126, 191, 233, 251, 255, 256},
		},
		{
			Total: 256,
			Dist:  []int{1, 5, 17, 45, 94, 157, 211, 241, 252, 255, 256},
		},
		{
			Total: 256,
			Dist:  []int{1, 5, 13, 33, 70, 125, 182, 223, 245, 253, 255, 256},
		},
		{
			Total: 256,
			Dist:  []int{1, 4, 11, 26, 54, 98, 151, 199, 232, 248, 254, 255, 256},
		},
		{
			Total: 256,
			Dist:  []int{1, 3, 9, 21, 42, 77, 124, 172, 212, 237, 249, 254, 255, 256},
		},
		{
			Total: 256,
			Dist:  []int{1, 2, 6, 16, 33, 60, 97, 144, 187, 220, 241, 250, 254, 255, 256},
		},
		{
			Total: 256,
			Dist:  []int{1, 2, 3, 11, 25, 47, 80, 120, 163, 201, 229, 245, 253, 254, 255, 256},
		},
		{
			Total: 256,
			Dist:  []int{1, 2, 3, 4, 17, 35, 62, 98, 139, 180, 214, 238, 252, 253, 254, 255, 256},
		},
	},
	{
		{
			Total: 256,
			Dist:  []int{127, 256},
		},
		{
			Total: 256,
			Dist:  []int{53, 202, 256},
		},
		{
			Total: 256,
			Dist:  []int{22, 127, 233, 256},
		},
		{
			Total: 256,
			Dist:  []int{11, 72, 183, 246, 256},
		},
		{
			Total: 256,
			Dist:  []int{6, 41, 127, 215, 251, 256},
		},
		{
			Total: 256,
			Dist:  []int{4, 24, 83, 170, 232, 253, 256},
		},
		{
			Total: 256,
			Dist:  []int{3, 16, 56, 127, 200, 241, 254, 256},
		},
		{
			Total: 256,
			Dist:  []int{3, 12, 39, 92, 162, 218, 246, 255, 256},
		},
		{
			Total: 256,
			Dist:  []int{3, 11, 30, 67, 124, 185, 229, 249, 255, 256},
		},
		{
			Total: 256,
			Dist:  []int{3, 10, 25, 53, 97, 151, 200, 233, 250, 255, 256},
		},
		{
			Total: 256,
			Dist:  []int{1, 8, 21, 43, 77, 123, 171, 209, 237, 251, 255, 256},
		},
		{
			Total: 256,
			Dist:  []int{1, 2, 13, 35, 62, 97, 139, 186, 219, 244, 254, 255, 256},
		},
		{
			Total: 256,
			Dist:  []int{1, 2, 8, 22, 48, 85, 128, 171, 208, 234, 248, 254, 255, 256},
		},
		{
			Total: 256,
			Dist:  []int{1, 2, 6, 16, 36, 67, 107, 149, 189, 220, 240, 250, 254, 255, 256},
		},
		{
			Total: 256,
			Dist:  []int{1, 2, 5, 13, 29, 55, 90, 128, 166, 201, 227, 243, 251, 254, 255, 256},
		},
		{
			Total: 256,
			Dist:  []int{1, 2, 4, 10, 22, 43, 73, 109, 147, 183, 213, 234, 246, 252, 254, 255, 256},
		},
	},
	{
		{
			Total: 256,
			Dist:  []int{127, 256},
		},
		{
			Total: 256,
			Dist:  []int{49, 206, 256},
		},
		{
			Total: 256,
			Dist:  []int{20, 127, 236, 256},
		},
		{
			Total: 256,
			Dist:  []int{11, 71, 184, 246, 256},
		},
		{
			Total: 256,
			Dist:  []int{7, 43, 127, 214, 250, 256},
		},
		{
			Total: 256,
			Dist:  []int{6, 30, 87, 169, 229, 252, 256},
		},
		{
			Total: 256,
			Dist:  []int{5, 23, 62, 126, 194, 236, 252, 256},
		},
		{
			Total: 256,
			Dist:  []int{6, 20, 49, 96, 157, 209, 239, 253, 256},
		},
		{
			Total: 256,
			Dist:  []int{1, 16, 39, 74, 125, 175, 215, 245, 255, 256},
		},
		{
			Total: 256,
			Dist:  []int{1, 2, 23, 55, 97, 149, 195, 236, 254, 255, 256},
		},
		{
			Total: 256,
			Dist:  []int{1, 7, 23, 50, 86, 128, 170, 206, 233, 249, 255, 256},
		},
		{
			Total: 256,
			Dist:  []int{1, 6, 18, 39, 70, 108, 148, 186, 217, 238, 250, 255, 256},
		},
		{
			Total: 256,
			Dist:  []int{1, 4, 13, 30, 56, 90, 128, 166, 200, 226, 243, 252, 255, 256},
		},
		{
			Total: 256,
			Dist:  []int{1, 4, 11, 25, 47, 76, 110, 146, 180, 209, 231, 245, 252, 255, 256},
		},
		{
			Total: 256,
			Dist:  []int{1, 3, 8, 19, 37, 62, 93, 128, 163, 194, 219, 237, 248, 253, 255, 256},
		},
		{
			Total: 256,
			Dist:  []int{1, 2, 6, 15, 30, 51, 79, 111, 145, 177, 205, 226, 241, 250, 254, 255, 256},
		},
	},
	{
		{
			Total: 256,
			Dist:  []int{128, 256},
		},
		{
			Total: 256,
			Dist:  []int{42, 214, 256},
		},
		{
			Total: 256,
			Dist:  []int{21, 128, 235, 256},
		},
		{
			Total: 256,
			Dist:  []int{12, 72, 184, 245, 256},
		},
		{
			Total: 256,
			Dist:  []int{8, 42, 128, 214, 249, 256},
		},
		{
			Total: 256,
			Dist:  []int{8, 31, 86, 176, 231, 251, 256},
		},
		{
			Total: 256,
			Dist:  []int{5, 20, 58, 130, 202, 238, 253, 256},
		},
		{
			Total: 256,
			Dist:  []int{6, 18, 45, 97, 174, 221, 241, 251, 256},
		},
		{
			Total: 256,
			Dist:  []int{6, 25, 53, 88, 128, 168, 203, 231, 250, 256},
		},
		{
			Total: 256,
			Dist:  []int{4, 18, 40, 71, 108, 148, 185, 216, 238, 252, 256},
		},
		{
			Total: 256,
			Dist:  []int{3, 13, 31, 57, 90, 128, 166, 199, 225, 243, 253, 256},
		},
		{
			Total: 256,
			Dist:  []int{2, 10, 23, 44, 73, 109, 147, 183, 212, 233, 246, 254, 256},
		},
		{
			Total: 256,
			Dist:  []int{1, 6, 16, 33, 58, 90, 128, 166, 198, 223, 240, 250, 255, 256},
		},
		{
			Total: 256,
			Dist:  []int{1, 5, 12, 25, 46, 75, 110, 146, 181, 210, 231, 244, 251, 255, 256},
		},
		{
			Total: 256,
			Dist:  []int{1, 3, 8, 18, 35, 60, 92, 128, 164, 196, 221, 238, 248, 253, 255, 256},
		},
		{
			Total: 256,
			Dist:  []int{1, 3, 7, 14, 27, 48, 76, 110, 146, 180, 208, 229, 242, 249, 253, 255, 256},
		},
	},
}

// EXC_LSB 激励最低有效位参数
var EXC_LSB = &comm.ICDFContext{
	Total: 256,
	Dist:  []int{136, 256},
}

// EXC_SIGN 激励符号参数
var EXC_SIGN = [][][]*comm.ICDFContext{
	{
		// Inactive 非活动状态
		{
			// Low offset 低偏移
			{
				Total: 256,
				Dist:  []int{2, 256},
			},
			{
				Total: 256,
				Dist:  []int{207, 256},
			},
			{
				Total: 256,
				Dist:  []int{189, 256},
			},
			{
				Total: 256,
				Dist:  []int{179, 256},
			},
			{
				Total: 256,
				Dist:  []int{174, 256},
			},
			{
				Total: 256,
				Dist:  []int{163, 256},
			},
			{
				Total: 256,
				Dist:  []int{157, 256},
			},
		},
		{
			// High offset 高偏移
			{
				Total: 256,
				Dist:  []int{58, 256},
			},
			{
				Total: 256,
				Dist:  []int{245, 256},
			},
			{
				Total: 256,
				Dist:  []int{238, 256},
			},
			{
				Total: 256,
				Dist:  []int{232, 256},
			},
			{
				Total: 256,
				Dist:  []int{225, 256},
			},
			{
				Total: 256,
				Dist:  []int{220, 256},
			},
			{
				Total: 256,
				Dist:  []int{211, 256},
			},
		},
	},
	{
		// Unvoiced 非浊音
		{
			// Low offset 低偏移
			{
				Total: 256,
				Dist:  []int{1, 256},
			},
			{
				Total: 256,
				Dist:  []int{210, 256},
			},
			{
				Total: 256,
				Dist:  []int{190, 256},
			},
			{
				Total: 256,
				Dist:  []int{178, 256},
			},
			{
				Total: 256,
				Dist:  []int{169, 256},
			},
			{
				Total: 256,
				Dist:  []int{162, 256},
			},
			{
				Total: 256,
				Dist:  []int{152, 256},
			},
		},
		{
			// High offset 高偏移
			{
				Total: 256,
				Dist:  []int{48, 256},
			},
			{
				Total: 256,
				Dist:  []int{242, 256},
			},
			{
				Total: 256,
				Dist:  []int{235, 256},
			},
			{
				Total: 256,
				Dist:  []int{224, 256},
			},
			{
				Total: 256,
				Dist:  []int{214, 256},
			},
			{
				Total: 256,
				Dist:  []int{205, 256},
			},
			{
				Total: 256,
				Dist:  []int{190, 256},
			},
		},
	},
	{
		// Voiced 浊音
		{
			// Low offset 低偏移
			{
				Total: 256,
				Dist:  []int{1, 256},
			},
			{
				Total: 256,
				Dist:  []int{162, 256},
			},
			{
				Total: 256,
				Dist:  []int{152, 256},
			},
			{
				Total: 256,
				Dist:  []int{147, 256},
			},
			{
				Total: 256,
				Dist:  []int{144, 256},
			},
			{
				Total: 256,
				Dist:  []int{141, 256},
			},
			{
				Total: 256,
				Dist:  []int{138, 256},
			},
		},
		{
			// High offset 高偏移
			{
				Total: 256,
				Dist:  []int{8, 256},
			},
			{
				Total: 256,
				Dist:  []int{203, 256},
			},
			{
				Total: 256,
				Dist:  []int{187, 256},
			},
			{
				Total: 256,
				Dist:  []int{176, 256},
			},
			{
				Total: 256,
				Dist:  []int{168, 256},
			},
			{
				Total: 256,
				Dist:  []int{161, 256},
			},
			{
				Total: 256,
				Dist:  []int{154, 256},
			},
		},
	},
}

// QUANT_OFFSET 量化偏移参数
var QUANT_OFFSET = [][]int32{
	{25, 60}, // Inactive or Unvoiced 非活动或非浊音
	{8, 25},  // Voiced 浊音
}

// SilkFrame 表示 SILK 音频帧
type SilkFrame struct {
	FrameType       FrameType
	LogGain         int
	Coded           bool
	PrevVoiced      bool
	Nlsfs           [16]int16
	Lpc             [16]float32
	InterpolatedLpc [16]float32
	Interpolated    bool
	InterpFactor4   bool
	PreviousLag     int32
	Output          []float32
	LpcHistory      []float32
}

// NewSilkFrame 创建新的 SilkFrame 实例
func NewSilkFrame() *SilkFrame {
	f := &SilkFrame{
		Output:     make([]float32, 2*LPC_HISTORY),
		LpcHistory: make([]float32, 2*LPC_HISTORY),
	}
	return f
}

// ParseSubframeGains 解析子帧增益
func (f *SilkFrame) ParseSubframeGains(rd *comm.RangeDecoder, coded bool) float32 {
	if coded {
		idx := f.FrameType.SignalTypeIndex()
		msb := int(rd.DecodeICDF(MSB_SUBFRAME_GAIN[idx]))
		lsb := int(rd.DecodeICDF(LSB_SUBFRAME_GAIN))
		f.LogGain = (msb<<3 | lsb)
		if f.LogGain < f.LogGain-16 {
			f.LogGain = f.LogGain - 16
		}
	} else {
		delta := int(rd.DecodeICDF(DELTA_SUBFRAME_GAIN))
		f.LogGain = delta*2 - 16
		if f.LogGain < f.LogGain+delta-4 {
			f.LogGain = f.LogGain + delta - 4
		}
		if f.LogGain < 0 {
			f.LogGain = 0
		} else if f.LogGain > 63 {
			f.LogGain = 63
		}
	}

	logGain := (f.LogGain * 0x1D1C71) >> 16
	logGain += 2090
	return float32(log2Lin(logGain)) / 65536.0
}

// ParseLpc 解析线性预测系数
func (f *SilkFrame) ParseLpc(rd *comm.RangeDecoder, interpolate bool, band Band) {
	idx := f.FrameType.VoicedIndex()
	lsfS1 := rd.DecodeICDF(band.Stage1[idx])

	// 获取相关参数
	mapCtx := band.Map[lsfS1]
	step := band.Step
	weightMap := band.PredWeight
	weightMapIndex := band.PredWeightIndex[lsfS1]
	weights := band.Weight[lsfS1]
	codebooks := band.Codebook[lsfS1]

	// 解码第二阶段线谱频率
	lsfsS2 := make([]int8, len(mapCtx))
	for i, icdf := range mapCtx {
		lsf := int8(rd.DecodeICDF(icdf)) - 4
		if lsf == -4 {
			lsf -= int8(rd.DecodeICDF(LSF_STAGE2_EXTENSION))
		} else if lsf == 4 {
			lsf += int8(rd.DecodeICDF(LSF_STAGE2_EXTENSION))
		}
		lsfsS2[i] = lsf
	}

	// 反量化步骤
	dequantStep := func(lsfS2 int16) int16 {
		fix := int32(0)
		if lsfS2 < 0 {
			fix = 102
		} else if lsfS2 > 0 {
			fix = -102
		}
		return int16((int32(lsfS2)*1024 + fix) * int32(step) >> 16)
	}

	// 计算残差
	var prev int16
	residuals := make([]int16, len(lsfsS2))
	for i := len(lsfsS2) - 1; i >= 0; i-- {
		ds := dequantStep(int16(lsfsS2[i]))
		res := ds
		if i != len(lsfsS2)-1 {
			weight := int32(weightMap[weightMapIndex[i]][i])
			res += int16((int32(prev) * weight) >> 8)
		}
		residuals[i] = res
		prev = res
	}

	// 计算归一化线谱频率
	nlsfs := make([]int16, len(residuals))
	for i := range residuals {
		r := residuals[len(residuals)-1-i]
		c := int32(codebooks[i])
		w := int32(weights[i])
		nlsf := (c << 7) + (int32(r)<<14)/w
		if nlsf < 0 {
			nlsf = 0
		} else if nlsf > 32767 {
			nlsf = 32767
		}
		nlsfs[i] = int16(nlsf)
	}

	// 稳定处理
	band.Stabilize(nlsfs)

	// 插值处理
	f.Interpolated = false
	f.InterpFactor4 = true
	if interpolate {
		weight := rd.DecodeICDF(LSF_INTERPOLATION_INDEX)
		if weight != 4 && f.Coded {
			f.Interpolated = true
			if weight != 0 {
				interpolatedNlsfs := make([]int16, len(nlsfs))
				for i := range nlsfs {
					prev := int32(f.Nlsfs[i])
					curr := int32(nlsfs[i])
					interpolatedNlsfs[i] = int16(prev + ((curr-prev)*int32(weight))>>2)
				}
				band.LsfToLpc(f.InterpolatedLpc[:], interpolatedNlsfs)
			} else {
				copy(f.InterpolatedLpc[:], f.Lpc[:])
			}
			f.InterpFactor4 = false
		}
	}

	// 更新LPC系数
	copy(f.Nlsfs[:], nlsfs)
	band.LsfToLpc(f.Lpc[:], nlsfs)
}

// ParsePitchLags 解析音高延迟
func (f *SilkFrame) ParsePitchLags(rd *comm.RangeDecoder, subframes []*SubFrame, absolute bool, pitch PitchLag) {
	// 解析绝对延迟
	parseAbsoluteLag := func(rd *comm.RangeDecoder) int32 {
		high := rd.DecodeICDF(PITCH_HIGH_PART)
		low := rd.DecodeICDF(pitch.LowPart())
		return int32(high)*int32(pitch.Scale()) + int32(low) + int32(pitch.MinLag())
	}

	var lag int32
	if !absolute {
		delta := rd.DecodeICDF(PITCH_DELTA)
		if delta != 0 {
			lag = f.PreviousLag + int32(delta) - 9
		} else {
			lag = parseAbsoluteLag(rd)
		}
	} else {
		lag = parseAbsoluteLag(rd)
	}

	f.PreviousLag = lag

	// 获取偏移量
	var offsets []int8
	if len(subframes) == 2 {
		idx := rd.DecodeICDF(pitch.Contour()[0])
		offsets = pitch.Offset()[0][idx]
	} else {
		idx := rd.DecodeICDF(pitch.Contour()[1])
		offsets = pitch.Offset()[1][idx]
	}

	// 设置每个子帧的音高延迟
	for i, sf := range subframes {
		offset := int32(offsets[i])
		sf.PitchLag = lag + offset
		if sf.PitchLag < int32(pitch.MinLag()) {
			sf.PitchLag = int32(pitch.MinLag())
		} else if sf.PitchLag > int32(pitch.MaxLag()) {
			sf.PitchLag = int32(pitch.MaxLag())
		}
	}
}

// ParseLtpFilterCoeff 解析长时预测滤波器系数
func (f *SilkFrame) ParseLtpFilterCoeff(rd *comm.RangeDecoder, subframes []*SubFrame) {
	idxPeriod := rd.DecodeICDF(LTP_PERIODICITY)

	for _, sf := range subframes {
		idxFilter := rd.DecodeICDF(LTP_FILTER[idxPeriod])
		filterTaps := LTP_TAPS[idxPeriod][idxFilter]
		for i, tap := range filterTaps {
			sf.LtpTaps[i] = float32(tap) / 128.0
		}
	}
}

// ParseExcitation 解析激励信号
func (f *SilkFrame) ParseExcitation(rd *comm.RangeDecoder, residuals []float32, longFrame bool, shell ShellBlock) {
	var shellBlocks int = 0
	if longFrame {
		shellBlocks = int(shell.ShellBlocks[1])
	} else {
		shellBlocks = int(shell.ShellBlocks[0])
	}

	pulseCount := make([]uint8, shellBlocks)
	lsbCount := make([]uint8, shellBlocks)
	excitation := make([]int32, shellBlocks*16)

	// 解码种子
	seed := uint32(rd.DecodeICDF(LCG_SEED))

	// 获取激励率
	voicedIndex := f.FrameType.VoicedIndex()
	rateLevel := rd.DecodeICDF(EXC_RATE[voicedIndex])

	// 解码脉冲计数
	for i := 0; i < int(shellBlocks); i++ {
		p := rd.DecodeICDF(PULSE_COUNT[rateLevel])
		if p == 17 {
			l := 0
			for l < 10 {
				l++
				p = rd.DecodeICDF(PULSE_COUNT[9])
				if p != 17 {
					break
				}
			}
			if l == 10 {
				p = rd.DecodeICDF(PULSE_COUNT[10])
			}
			lsbCount[i] = uint8(l)
		}
		pulseCount[i] = uint8(p)
	}

	// 解码激励位置
	for i := 0; i < shellBlocks; i++ {
		p := pulseCount[i]
		loc := excitation[i*16 : (i+1)*16]

		if p == 0 {
			for j := range loc {
				loc[j] = 0
			}
		} else {
			splitLoc := func(rd *comm.RangeDecoder, level int, avail int32) [2]int32 {
				if avail == 0 {
					return [2]int32{0, 0}
				}
				left := rd.DecodeICDF(PULSE_LOCATION[level][avail-1])
				return [2]int32{int32(left), int32(avail) - int32(left)}
			}

			dist := splitLoc(rd, 0, int32(p))
			for j := 0; j < 2; j++ {
				subLoc := loc[j*8 : (j+1)*8]
				subDist := splitLoc(rd, 1, dist[j])
				for k := 0; k < 2; k++ {
					subSubLoc := subLoc[k*4 : (k+1)*4]
					subSubDist := splitLoc(rd, 2, subDist[k])
					for l := 0; l < 2; l++ {
						subSubSubLoc := subSubLoc[l*2 : (l+1)*2]
						finalDist := splitLoc(rd, 3, subSubDist[l])
						subSubSubLoc[0] = finalDist[0]
						subSubSubLoc[1] = finalDist[1]
					}
				}
			}
		}
	}

	// 解码最低有效位
	for i := 0; i < shellBlocks; i++ {
		bits := lsbCount[i]
		loc := excitation[i*16 : (i+1)*16]
		for j := range loc {
			for b := 0; b < int(bits); b++ {
				loc[j] = (loc[j] << 1) | int32(rd.DecodeICDF(EXC_LSB))
			}
		}
	}

	// 解码符号
	for i := 0; i < shellBlocks; i++ {
		p := pulseCount[i]
		loc := excitation[i*16 : (i+1)*16]
		for j := range loc {
			if loc[j] != 0 {
				signalType := f.FrameType.SignalTypeIndex()
				qoffsetType := f.FrameType.QoffsetTypeIndex()
				pulse := int(p)
				if pulse > 6 {
					pulse = 6
				}

				sign := rd.DecodeICDF(EXC_SIGN[signalType][qoffsetType][pulse])
				if sign == 0 {
					loc[j] *= -1
				}
			}
		}
	}

	// 生成残差信号
	for i, l := range excitation {
		voiced := f.FrameType.VoicedIndex()
		qoffset := f.FrameType.QoffsetTypeIndex()
		offset := QUANT_OFFSET[voiced][qoffset]

		ex1 := l*256 + offset
		ex := ex1 - 20
		if l < 0 {
			ex = ex1 + 20
		}

		// 应用随机相位
		seed = seed*196314165 + 907633515
		if seed&0x80000000 != 0 {
			ex *= -1
		}
		seed += uint32(l)

		// 转换为浮点数
		residuals[i] = float32(ex) / 8388608.0
	}
} // Flush 重置帧状态
func (f *SilkFrame) Flush() {
	if f.Coded {
		f.LogGain = 0
		f.Coded = false
		f.PrevVoiced = false
		f.Nlsfs = [16]int16{}
		f.Lpc = [16]float32{}
		f.InterpolatedLpc = [16]float32{}
		f.Interpolated = false
		f.InterpFactor4 = false
		f.PreviousLag = 0

		// 清空并重置输出缓冲区
		f.Output = make([]float32, 2*LPC_HISTORY)
		f.LpcHistory = make([]float32, 2*LPC_HISTORY)
	}
}

// Parse 解析 SILK 帧
func (f *SilkFrame) Parse(rd *comm.RangeDecoder, info *SilkInfo, vad bool, first bool) error {
	// 解析帧类型
	if vad {
		switch rd.DecodeICDF(FRAME_TYPE_ACTIVE) {
		case 0:
			f.FrameType = FrameType{Active: true, Voiced: false, High: false}
		case 1:
			f.FrameType = FrameType{Active: true, Voiced: false, High: true}
		case 2:
			f.FrameType = FrameType{Active: true, Voiced: true, High: false}
		case 3:
			f.FrameType = FrameType{Active: true, Voiced: true, High: true}
		}
	} else {
		if rd.DecodeICDF(FRAME_TYPE_INACTIVE) == 0 {
			f.FrameType = FrameType{Active: false, Voiced: false, High: false}
		} else {
			f.FrameType = FrameType{Active: false, Voiced: false, High: true}
		}
	}

	// 创建子帧和残差缓冲区
	subframes := make([]*SubFrame, info.Subframes)
	for i := range subframes {
		subframes[i] = &SubFrame{}
	}
	residuals := make([]float32, LPC_HISTORY+RES_HISTORY)

	// 解析子帧增益
	for i, sf := range subframes {
		coded := i == 0 && (first || !f.Coded)
		sf.Gain = f.ParseSubframeGains(rd, coded)
	}

	// 确定是否为长帧
	longFrame := info.Subframes == 4

	// 解析 LPC 系数
	var order int
	if info.Bandwidth > BandwidthMedium {
		f.ParseLpc(rd, longFrame, WB)
		order = WB.Order
	} else {
		f.ParseLpc(rd, longFrame, NB_MB)
		order = NB_MB.Order
	}

	// 解析音高延迟和 LTP 系数
	if f.FrameType.Voiced {
		absolute := first || !f.PrevVoiced
		switch info.Bandwidth {
		case BandwidthNarrow:
			f.ParsePitchLags(rd, subframes, absolute, NB_PitchLag)
		case BandwidthMedium:
			f.ParsePitchLags(rd, subframes, absolute, MB_PitchLag)
		default:
			f.ParsePitchLags(rd, subframes, absolute, WB_PitchLag)
		}
		f.ParseLtpFilterCoeff(rd, subframes)
	}

	// 计算 LTP 缩放因子
	ltpScale := float32(15565.0 / 16384.0)
	if f.FrameType.Voiced && first {
		idx := rd.DecodeICDF(LTP_SCALE_INDEX)
		ltpScale = float32(LTP_SCALE[idx]) / 16384.0
	}

	// 解析激励信号
	switch info.Bandwidth {
	case BandwidthNarrow:
		f.ParseExcitation(rd, residuals[RES_HISTORY:], longFrame, NBShellBlock)
	case BandwidthMedium:
		f.ParseExcitation(rd, residuals[RES_HISTORY:], longFrame, MBShellBlock)
	default:
		f.ParseExcitation(rd, residuals[RES_HISTORY:], longFrame, WBShellBlock)
	}

	// 处理每个子帧
	for i, sf := range subframes {
		// 选择 LPC 系数
		var lpcCoeff []float32
		if i < 2 && f.Interpolated {
			lpcCoeff = f.InterpolatedLpc[:order]
		} else {
			lpcCoeff = f.Lpc[:order]
		}

		// 浊音帧处理
		if f.FrameType.Voiced {
			before := int(sf.PitchLag) + LTP_ORDER/2
			end := 0
			scale := float32(1.0)

			if i < 2 || f.InterpFactor4 {
				end = i * info.SfSize
				scale = ltpScale
			} else {
				end = (i - 2) * info.SfSize
			}

			// 重新白化残差
			if before > end {
				start := LPC_HISTORY + i*info.SfSize - before
				stop := LPC_HISTORY + i*info.SfSize - end
				startRes := RES_HISTORY + i*info.SfSize - before
				//stopRes := RES_HISTORY + i*info.SfSize - end

				for j := start; j < stop; j++ {
					// 获取历史窗口
					history := f.Output[j-order : j]

					// 计算预测值
					sum := f.Output[j]
					for k, coeff := range lpcCoeff {
						sum -= coeff * history[len(history)-1-k]
					}

					// 更新残差
					residuals[startRes+j-start] = clamp(sum, -1, 1) * scale / sf.Gain
				}
			}

			// 重新缩放残差
			if end != 0 {
				start := RES_HISTORY + i*info.SfSize - end
				stop := RES_HISTORY + i*info.SfSize
				rescale := subframes[i-1].Gain / sf.Gain

				for j := start; j < stop; j++ {
					residuals[j] *= rescale
				}
			}

			// 应用长时预测
			start := RES_HISTORY + i*info.SfSize
			stop := start + info.SfSize

			for j := start; j < stop; j++ {
				sum := residuals[j]
				for k := 0; k < LTP_ORDER; k++ {
					idx := j - int(sf.PitchLag) + LTP_ORDER/2 - k
					sum += sf.LtpTaps[k] * residuals[idx]
				}
				residuals[j] = sum
			}
		}

		// 合成音频
		startLpc := LPC_HISTORY + i*info.SfSize
		//stopLpc := LPC_HISTORY + (i+1)*info.SfSize
		startRes := RES_HISTORY + i*info.SfSize

		for j := 0; j < info.SfSize; j++ {
			// 计算激励
			sum := residuals[startRes+j] * sf.Gain

			// 应用 LPC 合成
			for k := 0; k < order; k++ {
				idx := startLpc + j - k - 1
				if idx < 0 {
					idx += len(f.LpcHistory)
				}
				sum += lpcCoeff[k] * f.LpcHistory[idx]
			}

			// 更新历史并存储输出
			f.LpcHistory[startLpc+j] = sum
			f.Output[startLpc+j] = clamp(sum, -1, 1)
		}
	}

	// 更新帧状态
	f.PrevVoiced = f.FrameType.Voiced
	f.Coded = true

	// 移动历史缓冲区
	for i := 0; i < LPC_HISTORY; i++ {
		f.LpcHistory[i] = f.LpcHistory[i+info.FSize]
		f.Output[i] = f.Output[i+info.FSize]
	}

	return nil
}

// clamp 限制值在[min, max]范围内
func clamp(value, min, max float32) float32 {
	if value < min {
		return min
	}
	if value > max {
		return max
	}
	return value
}

// NewSilk 创建新的 SILK 解码器
func NewSilk(stereoOut bool) *Silk {
	return &Silk{
		StereoOut: stereoOut,
		MidFrame:  NewSilkFrame(),
		SideFrame: NewSilkFrame(),
	}
}

// Flush 重置解码器状态
func (s *Silk) Flush() {
	s.MidFrame.Flush()
	s.SideFrame.Flush()

	s.Info.Prev0 = 0.0
	s.Info.Prev1 = 0.0
}

// Setup 配置解码器参数
func (s *Silk) Setup(pkt *comm.Packet) {
	switch pkt.FrameDuration {
	case comm.FrameDurationMedium:
		s.Frames = 1
		s.Info.Subframes = 2
	case comm.FrameDurationStandard:
		s.Frames = 1
		s.Info.Subframes = 4
	case comm.FrameDurationLong:
		s.Frames = 2
		s.Info.Subframes = 4
	case comm.FrameDurationVeryLong:
		s.Frames = 3
		s.Info.Subframes = 4
	}
	s.Stereo = pkt.Stereo
	s.Info.Bandwidth = Bandwidth(pkt.Bandwidth)
	if s.Info.Bandwidth > BandwidthWide {
		s.Info.Bandwidth = BandwidthWide
	}

	switch s.Info.Bandwidth {
	case BandwidthNarrow:
		s.Info.SfSize = 40
	case BandwidthMedium:
		s.Info.SfSize = 60
	case BandwidthWide:
		s.Info.SfSize = 80
	}
	s.Info.FSize = s.Info.SfSize * s.Info.Subframes

	// 重置输出缓冲区
	s.LeftOutbuf = make([]float32, s.Info.FSize*s.Frames)
	s.RightOutbuf = make([]float32, s.Info.FSize*s.Frames)
}

// ParseStereoWeight 解析立体声权重
func (s *Silk) ParseStereoWeight(rd *comm.RangeDecoder, vad bool) bool {
	wQ13 := []int16{
		-13732, -10050, -8266, -7526, -6500, -5000, -2950, -820,
		820, 2950, 5000, 6500, 7526, 8266, 10050, 13732,
	}

	n := rd.DecodeICDF(STAGE1)
	i0 := rd.DecodeICDF(STAGE2) + 3*(n/5)
	i1 := rd.DecodeICDF(STAGE3)*2 + 1
	i2 := rd.DecodeICDF(STAGE2) + 3*(n%5)
	i3 := rd.DecodeICDF(STAGE3)*2 + 1

	weight := func(idx, scale int) int16 {
		w := wQ13[idx]
		w1 := wQ13[idx+1]
		return w + int16((int32(w1-w)*6554)>>16*int32(scale))
	}

	w0 := weight(int(i0), int(i1))
	w1 := weight(int(i2), int(i3))

	s.Info.Weight0 = float32(w0-w1) / 8192.0
	s.Info.Weight1 = float32(w1) / 8192.0

	if vad {
		return false
	}
	return rd.DecodeICDF(MID_ONLY) != 0
}

// UnmixMs 解混中-侧声道
func (s *Silk) UnmixMs(start, end int) {
	inStart := LPC_HISTORY - s.Info.FSize
	inRange := inStart + s.Info.FSize
	w0 := s.Info.Weight0
	w1 := s.Info.Weight1
	w0p := s.Info.Prev0
	w1p := s.Info.Prev1

	n1 := 128
	switch s.Info.Bandwidth {
	case BandwidthNarrow:
		n1 = 64
	case BandwidthMedium:
		n1 = 96
	}

	w0d := (w0 - w0p) / float32(n1)
	w1d := (w1 - w1p) / float32(n1)

	mid := s.MidFrame.Output[inStart-2 : inRange]
	side := s.SideFrame.Output[inStart-1 : inRange-1]

	for i := 0; i < s.Info.FSize; i++ {
		idx := start + i
		var interp0, interp1 float32
		if i < n1 {
			interp0 = w0p + float32(i)*w0d
			interp1 = w1p + float32(i)*w1d
		} else {
			interp0 = w0
			interp1 = w1
		}

		// 计算预测值
		var p0 float32
		if i < s.Info.FSize-2 {
			p0 = 0.25 * (mid[i] + 2*mid[i+1] + mid[i+2])
		} else {
			p0 = mid[i+1]
		}

		si0 := side[i] + interp0*p0
		r := (1.0+interp1)*mid[i+1] + si0
		l := (1.0-interp1)*mid[i+1] - si0

		// 限制在[-1,1]范围内
		if r > 1.0 {
			r = 1.0
		} else if r < -1.0 {
			r = -1.0
		}
		if l > 1.0 {
			l = 1.0
		} else if l < -1.0 {
			l = -1.0
		}

		s.RightOutbuf[idx] = r
		s.LeftOutbuf[idx] = l
	}

	s.Info.Prev0 = s.Info.Weight0
	s.Info.Prev1 = s.Info.Weight1
}

// Decode 解码 SILK 音频帧
func (s *Silk) Decode(rd *comm.RangeDecoder) (int, error) {
	midVad := make([]bool, s.Frames)
	sideVad := make([]bool, s.Frames)

	// 解析语音活动检测
	lp := func(rd *comm.RangeDecoder, vad []bool) error {
		for i := range vad {
			vad[i] = rd.DecodeLogP(1)
		}
		if rd.DecodeLogP(1) {
			return errors.New("unsupported LBRR frames")
		}
		return nil
	}

	if err := lp(rd, midVad[:s.Frames]); err != nil {
		return 0, err
	}

	if s.Stereo {
		if err := lp(rd, sideVad[:s.Frames]); err != nil {
			return 0, err
		}
	}

	for i := 0; i < s.Frames; i++ {
		first := i == 0
		midOnly := false
		if s.Stereo {
			midOnly = s.ParseStereoWeight(rd, sideVad[i])
		}

		// 解析中声道
		if err := s.MidFrame.Parse(rd, &s.Info, midVad[i], first); err != nil {
			return 0, err
		}

		// 解析侧声道
		if s.Stereo && !midOnly {
			if err := s.SideFrame.Parse(rd, &s.Info, sideVad[i], first); err != nil {
				return 0, err
			}
		} else if s.Stereo {
			s.SideFrame.Flush()
		}

		// 生成输出
		outStart := i * s.Info.FSize
		outEnd := (i + 1) * s.Info.FSize
		//outRange := outStart:outEnd

		if s.Stereo && s.StereoOut {
			s.UnmixMs(outStart, outEnd)
		} else {
			inStart := LPC_HISTORY - s.Info.FSize - 2
			inEnd := inStart + s.Info.FSize
			midOutput := s.MidFrame.Output[inStart:inEnd]

			if s.StereoOut {
				copy(s.LeftOutbuf[outStart:outEnd], midOutput)
			}
			copy(s.RightOutbuf[outStart:outEnd], midOutput)
		}
	}

	return s.Info.FSize * s.Frames, nil
}
