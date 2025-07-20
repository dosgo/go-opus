package celt

import (
	"fmt"
	"math"

	"github.com/dosgo/goOpus/comm"
)

const (
	SHORT_BLOCKSIZE        = 120
	OVERLAP                = SHORT_BLOCKSIZE
	MAX_LOG_BLOCKS         = 3
	MAX_FRAME_SIZE         = SHORT_BLOCKSIZE * (1 << MAX_LOG_BLOCKS)
	MAX_BANDS              = 21
	MIN_PERIOD             = 15
	SPREAD_NONE            = 0
	SPREAD_LIGHT           = 1
	SPREAD_NORMAL          = 2
	SPREAD_AGGRESSIVE      = 3
	FRAC_1_SQRT_2          = 0.7071067811865476
	MAX_FINE_BITS          = 8
	QTHETA_OFFSET          = 4
	QTHETA_OFFSET_TWOPHASE = 16
	BITRES                 = 16
)

const MAX_FRAMES = 2
const ALLOC_STEPS = 6
const CELT_VECTOR = 11

/*
	type PostFilter struct {
		Period    int
		PeriodNew int
		PeriodOld int
		Gains     [3]float32
		GainsNew  [3]float32
		GainsOld  [3]float32
	}
*/
type PostFilter struct {
	period     int
	period_new int
	period_old int
	gains      [3]float32
	gains_new  [3]float32
	gains_old  [3]float32
}

/*
	type CeltFrame struct {
		PF           PostFilter
		Energy       [MAX_BANDS]float32
		PrevEnergy   [MAX_BANDS]float32
		CollapseMask [MAX_BANDS]uint8
		Buf          [2048]float32
		DeemphCoeff  float32
	}
*/
type CeltFrame struct {
	pf             PostFilter
	energy         [MAX_BANDS]float32
	prev_energy    [MAX_BANDS]float32
	collapse_masks [MAX_BANDS]uint8
	buf            [2048]float32
	deemph_coeff   float32
}

/*
type Celt struct {
	Stereo          bool
	StereoPkt       bool
	Bits            int
	Lm              int
	BandStart       int
	BandEnd         int
	Frames          [2]CeltFrame
	Spread          int
	FineBits        [MAX_BANDS]int32
	FinePriority    [MAX_BANDS]bool
	Pulses          [MAX_BANDS]int32
	TfChange        [MAX_BANDS]int8
	AnticollapseBit int
	Blocks          int
	Blocksize       int
	IntensityStereo int
	DualStereo      bool
	Remaining       int32
	Remaining2      int32
	Codedband       int
	Scratch         [22 * 8]float32
	Seed            uint32
}*/
// Celt 结构体
type Celt struct {
	stereo     bool
	stereo_pkt bool
	bits       int // usize -> int (平台相关)
	lm         int // usize -> int
	band_start int // Range<usize> 起始
	band_end   int // Range<usize> 结束
	frames     [2]CeltFrame
	spread     int // usize -> int

	fine_bits     [MAX_BANDS]int32
	fine_priority [MAX_BANDS]bool
	pulses        [MAX_BANDS]int32
	tf_change     [MAX_BANDS]int8

	anticollapse_bit int // usize -> int
	blocks           int // usize -> int
	blocksize        int // usize -> int

	intensity_stereo int // usize -> int
	dual_stereo      bool

	remaining  int32
	remaining2 int32
	codedband  int // usize -> int

	scratch [22 * 8]float32
	seed    uint32
}

var POSTFILTER_TAPS = [][3]float32{
	{0.3066406250, 0.2170410156, 0.1296386719},
	{0.4638671875, 0.2680664062, 0.0},
	{0.7998046875, 0.1000976562, 0.0},
}
var TAPSET = &comm.ICDFContext{
	Total: 4,
	Dist:  []int{2, 3, 4},
}
var ALPHA_COEF = []float32{
	29440.0 / 32768.0,
	26112.0 / 32768.0,
	21248.0 / 32768.0,
	16384.0 / 32768.0,
}

var BETA_COEF = []float32{
	1.0 - 30147.0/32768.0,
	1.0 - 22282.0/32768.0,
	1.0 - 12124.0/32768.0,
	1.0 - 6554.0/32768.0,
}

var COARSE_ENERGY_INTRA = [][]uint8{
	// 120-samples
	{
		24, 179, 48, 138, 54, 135, 54, 132, 53, 134, 56, 133, 55, 132, 55, 132,
		61, 114, 70, 96, 74, 88, 75, 88, 87, 74, 89, 66, 91, 67, 100, 59,
		108, 50, 120, 40, 122, 37, 97, 43, 78, 50,
	},
	// 240-samples
	{
		23, 178, 54, 115, 63, 102, 66, 98, 69, 99, 74, 89, 71, 91, 73, 91,
		78, 89, 86, 80, 92, 66, 93, 64, 102, 59, 103, 60, 104, 60, 117, 52,
		123, 44, 138, 35, 133, 31, 97, 38, 77, 45,
	},
	// 480-samples
	{
		21, 178, 59, 110, 71, 86, 75, 85, 84, 83, 91, 66, 88, 73, 87, 72,
		92, 75, 98, 72, 105, 58, 107, 54, 115, 52, 114, 55, 112, 56, 129, 51,
		132, 40, 150, 33, 140, 29, 98, 35, 77, 42,
	},
	// 960-samples
	{
		22, 178, 63, 114, 74, 82, 84, 83, 92, 82, 103, 62, 96, 72, 96, 67,
		101, 73, 107, 72, 113, 55, 118, 52, 125, 52, 118, 52, 117, 55, 135, 49,
		137, 39, 157, 32, 145, 29, 97, 33, 77, 40,
	},
}

var COARSE_ENERGY_INTER = [][]uint8{
	// 120-samples
	{
		72, 127, 65, 129, 66, 128, 65, 128, 64, 128, 62, 128, 64, 128, 64, 128,
		92, 78, 92, 79, 92, 78, 90, 79, 116, 41, 115, 40, 114, 40, 132, 26,
		132, 26, 145, 17, 161, 12, 176, 10, 177, 11,
	},
	// 240-samples
	{
		83, 78, 84, 81, 88, 75, 86, 74, 87, 71, 90, 73, 93, 74, 93, 74,
		109, 40, 114, 36, 117, 34, 117, 34, 143, 17, 145, 18, 146, 19, 162, 12,
		165, 10, 178, 7, 189, 6, 190, 8, 177, 9,
	},
	// 480-samples
	{
		61, 90, 93, 60, 105, 42, 107, 41, 110, 45, 116, 38, 113, 38, 112, 38,
		124, 26, 132, 27, 136, 19, 140, 20, 155, 14, 159, 16, 158, 18, 170, 13,
		177, 10, 187, 8, 192, 6, 175, 9, 159, 10,
	},
	// 960-samples
	{
		42, 121, 96, 66, 108, 43, 111, 40, 117, 44, 123, 32, 120, 36, 119, 33,
		127, 33, 134, 34, 139, 21, 147, 23, 152, 20, 158, 25, 154, 26, 166, 21,
		173, 16, 184, 13, 184, 10, 150, 13, 139, 15,
	},
}

var STATIC_CAPS = [][][]uint8{
	// 120-sample
	{
		{
			224, 224, 224, 224, 224, 224, 224, 224, 160, 160, 160, 160, 185, 185, 185, 178, 178,
			168, 134, 61, 37,
		},
		{
			224, 224, 224, 224, 224, 224, 224, 224, 240, 240, 240, 240, 207, 207, 207, 198, 198,
			183, 144, 66, 40,
		},
	},
	// 240-sample
	{
		{
			160, 160, 160, 160, 160, 160, 160, 160, 185, 185, 185, 185, 193, 193, 193, 183, 183,
			172, 138, 64, 38,
		},
		{
			240, 240, 240, 240, 240, 240, 240, 240, 207, 207, 207, 207, 204, 204, 204, 193, 193,
			180, 143, 66, 40,
		},
	},
	// 480-sample
	{
		{
			185, 185, 185, 185, 185, 185, 185, 185, 193, 193, 193, 193, 193, 193, 193, 183, 183,
			172, 138, 65, 39,
		},
		{
			207, 207, 207, 207, 207, 207, 207, 207, 204, 204, 204, 204, 201, 201, 201, 188, 188,
			176, 141, 66, 40,
		},
	},
	// 960-sample
	{
		{
			193, 193, 193, 193, 193, 193, 193, 193, 193, 193, 193, 193, 194, 194, 194, 184, 184,
			173, 139, 65, 39,
		},
		{
			204, 204, 204, 204, 204, 204, 204, 204, 201, 201, 201, 201, 198, 198, 198, 187, 187,
			175, 140, 66, 40,
		},
	},
}
var FREQ_RANGE = []uint8{
	1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 4, 4, 4, 6, 6, 8, 12, 18, 22,
}

var MODEL_ENERGY_SMALL = &comm.ICDFContext{
	Total: 4,
	Dist:  []int{2, 3, 4},
}

var TF_SELECT = [][][][2]int8{
	{
		{{0, -1}, {0, -1}},
		{{0, -1}, {0, -1}},
	},
	{
		{{0, -1}, {0, -2}},
		{{1, 0}, {1, -1}},
	},
	{
		{{0, -2}, {0, -3}},
		{{2, 0}, {1, -1}},
	},
	{
		{{0, -2}, {0, -3}},
		{{3, 0}, {1, -1}},
	},
}

var MODEL_SPREAD = &comm.ICDFContext{
	Total: 32,
	Dist:  []int{7, 9, 30, 32},
}

var ALLOC_TRIM = &comm.ICDFContext{
	Total: 128,
	Dist:  []int{2, 4, 9, 19, 41, 87, 109, 119, 124, 126, 128},
}

var LOG2_FRAC = []uint8{
	0, 8, 13, 16, 19, 21, 23, 24, 26, 27, 28, 29, 30, 31, 32, 32, 33, 34, 34, 35, 36, 36, 37, 37,
}

var STATIC_ALLOC = [][MAX_BANDS]uint8{
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{90, 80, 75, 69, 63, 56, 49, 40, 34, 29, 20, 18, 10, 0, 0, 0, 0, 0, 0, 0, 0},
	{110, 100, 90, 84, 78, 71, 65, 58, 51, 45, 39, 32, 26, 20, 12, 0, 0, 0, 0, 0, 0},
	{118, 110, 103, 93, 86, 80, 75, 70, 65, 59, 53, 47, 40, 31, 23, 15, 4, 0, 0, 0, 0},
	{126, 119, 112, 104, 95, 89, 83, 78, 72, 66, 60, 54, 47, 39, 32, 25, 17, 12, 1, 0, 0},
	{134, 127, 120, 114, 103, 97, 91, 85, 78, 72, 76, 70, 64, 57, 51, 45, 39, 33, 26, 15, 1},
	{144, 137, 130, 124, 113, 107, 101, 95, 88, 82, 76, 70, 64, 57, 51, 45, 39, 33, 26, 15, 1},
	{152, 145, 138, 132, 123, 117, 111, 105, 98, 92, 86, 80, 74, 67, 61, 55, 49, 43, 36, 20, 1},
	{162, 155, 148, 142, 133, 127, 121, 115, 108, 102, 96, 90, 84, 77, 71, 65, 59, 53, 46, 30, 1},
	{172, 165, 158, 152, 143, 137, 131, 125, 118, 112, 106, 100, 94, 87, 81, 75, 69, 63, 56, 45, 20},
	{200, 200, 200, 200, 200, 200, 200, 200, 198, 193, 188, 183, 178, 173, 168, 163, 158, 153, 148, 129, 104},
}

var FREQ_BANDS = []uint8{
	0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 20, 24, 28, 34, 40, 48, 60, 78, 100,
}

var LOG_FREQ_RANGE = []uint8{
	0, 0, 0, 0, 0, 0, 0, 0, 8, 8, 8, 8, 16, 16, 16, 21, 21, 24, 29, 34, 36,
}

var BIT_INTERLEAVE = []uint8{
	0, 1, 1, 1, 2, 3, 3, 3, 2, 3, 3, 3, 2, 3, 3, 3,
}

var PVQ_U = []uint32{
	// PVQ_U数据太大，这里只展示部分
	1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51,
	13, 25, 41, 61, 85, 113, 145, 181, 221, 265, 313, 365, 421, 481, 545, 613, 685, 761, 841, 925,
	63, 129, 231, 377, 575, 833, 1159, 1561, 2047, 2625, 3303, 4089, 4991, 6017, 7175, 8473, 9919,
	1683, 3653, 7183, 13073, 22363, 36365, 56695, 85305, 124515, 177045, 246047, 335137, 448427,
	8989, 19825, 40081, 75517, 134245, 227305, 369305, 579125, 880685, 1303777, 1884961, 2653649,
	48639, 108545, 224143, 433905, 795455, 1392065, 2340495, 3800305, 5984767, 9173505, 13726991,
	265729, 598417, 1256465, 2485825, 4673345, 8405905, 14546705, 24331777, 39490049, 62390545,
	1462563, 3317445, 7059735, 14218905, 27298155, 50250765, 89129247, 152951073, 254831667,
	8097453, 18474633, 39753273, 81270333, 158819253, 298199265, 540279585, 948062325,
	45046719, 103274625, 224298231, 464387817, 921406335, 1759885185, 3248227095,
	251595969, 579168825, 1267854873, 2653649025,
	1409933619,
}

var PVQ_U_ROW = []uint32{0, 176, 351, 525, 698, 870, 1041, 1131, 1178, 1207, 1226, 1240, 1248, 1254, 1257}

func pvqURow(rowIndex int) []uint32 {
	start := PVQ_U_ROW[rowIndex]
	return PVQ_U[start:]
}

var CACHE_BITS = []uint8{
	40, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
	7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 40, 15, 23, 28, 31, 34, 36, 38, 39, 41, 42, 43, 44, 45, 46, 47,
	47, 49, 50, 51, 52, 53, 54, 55, 55, 57, 58, 59, 60, 61, 62, 63, 63, 65, 66, 67, 68, 69, 70, 71,
	71, 40, 20, 33, 41, 48, 53, 57, 61, 64, 66, 69, 71, 73, 75, 76, 78, 80, 82, 85, 87, 89, 91, 92,
	94, 96, 98, 101, 103, 105, 107, 108, 110, 112, 114, 117, 119, 121, 123, 124, 126, 128, 40, 23,
	39, 51, 60, 67, 73, 79, 83, 87, 91, 94, 97, 100, 102, 105, 107, 111, 115, 118, 121, 124, 126,
	129, 131, 135, 139, 142, 145, 148, 150, 153, 155, 159, 163, 166, 169, 172, 174, 177, 179, 35,
	28, 49, 65, 78, 89, 99, 107, 114, 120, 126, 132, 136, 141, 145, 149, 153, 159, 165, 171, 176,
	180, 185, 189, 192, 199, 205, 211, 216, 220, 225, 229, 232, 239, 245, 251, 21, 33, 58, 79, 97,
	112, 125, 137, 148, 157, 166, 174, 182, 189, 195, 201, 207, 217, 227, 235, 243, 251, 17, 35,
	63, 86, 106, 123, 139, 152, 165, 177, 187, 197, 206, 214, 222, 230, 237, 250, 25, 31, 55, 75,
	91, 105, 117, 128, 138, 146, 154, 161, 168, 174, 180, 185, 190, 200, 208, 215, 222, 229, 235,
	240, 245, 255, 16, 36, 65, 89, 110, 128, 144, 159, 173, 185, 196, 207, 217, 226, 234, 242, 250,
	11, 41, 74, 103, 128, 151, 172, 191, 209, 225, 241, 255, 9, 43, 79, 110, 138, 163, 186, 207,
	227, 246, 12, 39, 71, 99, 123, 144, 164, 182, 198, 214, 228, 241, 253, 9, 44, 81, 113, 142,
	168, 192, 214, 235, 255, 7, 49, 90, 127, 160, 191, 220, 247, 6, 51, 95, 134, 170, 203, 234, 7,
	47, 87, 123, 155, 184, 212, 237, 6, 52, 97, 137, 174, 208, 240, 5, 57, 106, 151, 192, 231, 5,
	59, 111, 158, 202, 243, 5, 55, 103, 147, 187, 224, 5, 60, 113, 161, 206, 248, 4, 65, 122, 175,
	224, 4, 67, 127, 182, 234,
}

var CACHE_INDEX = []int16{
	-1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 41, 41, 41, 82, 82, 123, 164, 200, 222, 0, 0, 0, 0,
	0, 0, 0, 0, 41, 41, 41, 41, 123, 123, 123, 164, 164, 240, 266, 283, 295, 41, 41, 41, 41, 41,
	41, 41, 41, 123, 123, 123, 123, 240, 240, 240, 266, 266, 305, 318, 328, 336, 123, 123, 123,
	123, 123, 123, 123, 123, 240, 240, 240, 240, 305, 305, 305, 318, 318, 343, 351, 358, 364, 240,
	240, 240, 240, 240, 240, 240, 240, 305, 305, 305, 305, 343, 343, 343, 351, 351, 370, 376, 382,
	387,
}

var QN_EXP2 = []uint16{16384, 17866, 19483, 21247, 23170, 25267, 27554, 30048}

var BIT_DEINTERLEAVE = []uint8{
	0x00, 0x03, 0x0C, 0x0F, 0x30, 0x33, 0x3C, 0x3F, 0xC0, 0xC3, 0xCC, 0xCF, 0xF0, 0xF3, 0xFC, 0xFF,
}

const invSqrt2 = float32(0.7071067811865476) // 1 / math.Sqrt(2)

func haar1(buf []float32, n0, stride int) {
	for i := 0; i < n0/2; i++ {
		offset := i * 2 * stride
		l0 := buf[offset : offset+stride]
		l1 := buf[offset+stride : offset+2*stride]

		for j := 0; j < stride; j++ {
			e0 := l0[j]
			e1 := l1[j]
			v0 := (e0 + e1) * invSqrt2
			v1 := (e0 - e1) * invSqrt2
			l0[j] = v0
			l1[j] = v1
		}
	}
}

var HADAMARD_ORDERY = []int{
	1, 0, 3, 0, 2, 1, 7, 0, 4, 3, 6, 1, 5, 2, 15, 0, 8, 7, 12, 3, 11, 4, 14, 1, 9, 6, 13, 2, 10, 5,
}

func interleaveHadamard(scratch, buf []float32, n0, stride int, hadamard bool) {
	size := n0 * stride

	if hadamard {
		shuffle := HADAMARD_ORDERY[:stride] // 注意：原 Rust 代码使用 stride-2，但这里调整索引
		for i := 0; i < stride; i++ {
			for j := 0; j < n0; j++ {
				scratch[j*stride+i] = buf[shuffle[i]*n0+j]
			}
		}
	} else {
		for i := 0; i < stride; i++ {
			for j := 0; j < n0; j++ {
				scratch[j*stride+i] = buf[i*n0+j]
			}
		}
	}

	fmt.Println("interleave")
	for _, v := range buf[:size] {
		fmt.Printf("  %#.10g\n", v)
	}

	copy(buf[:size], scratch[:size])
}

func deinterleaveHadamard(scratch, buf []float32, n0, stride int, hadamard bool) {
	size := n0 * stride

	fmt.Println("before deinterleave")
	for _, v := range buf[:size] {
		fmt.Printf("  %#.10g\n", v)
	}

	if hadamard {
		shuffle := HADAMARD_ORDERY[:stride]
		for i := 0; i < stride; i++ {
			for j := 0; j < n0; j++ {
				scratch[shuffle[i]*n0+j] = buf[j*stride+i]
			}
		}
	} else {
		for i := 0; i < stride; i++ {
			for j := 0; j < n0; j++ {
				scratch[i*n0+j] = buf[j*stride+i]
			}
		}
	}

	fmt.Println("deinterleave")
	for _, v := range scratch[:size] {
		fmt.Printf("  %#.10g\n", v)
	}

	copy(buf[:size], scratch[:size])
}

// 辅助函数
func update(k0, k uint32, s int32, norm *uint32) int32 {
	fmt.Printf("%d - %d\n", k0, k)
	d := int64(k0) - int64(k)
	val := int32(d) + s ^ s
	*norm += uint32(val * val)
	return val
}

func cwrsi(n, k, i uint32, y []int32) uint32 {
	norm := uint32(0)
	yPtr := 0

	for n > 2 {
		fmt.Printf("k %d n %d i %d\n", k, n, i)
		if k >= n {
			row := pvqURow(int(n))
			p := uint32(row[k+1])
			fmt.Printf("pulse %d\n", p)

			s := int32(0)
			if i >= p {
				i -= p
				s = -1
			}

			k0 := k
			q := uint32(row[n])
			var pVal uint32

			if q > i {
				k = n
				for {
					k--
					rowK := pvqURow(int(k))
					pVal = uint32(rowK[n])
					fmt.Printf("pulse %d\n", pVal)
					if i >= pVal {
						break
					}
				}
			} else {
				for {
					rowK := pvqURow(int(k))
					pVal = uint32(rowK[k])
					if i >= pVal {
						break
					}
					k--
				}
			}

			i -= pVal
			fmt.Printf("-- i %d p %d\n", i, pVal)
			y[yPtr] = update(k0, k, s, &norm)
			yPtr++
		} else {
			row := pvqURow(int(k))
			p := uint32(row[n])
			q := uint32(row[k+1])
			fmt.Printf("i %d p %d q %d\n", i, p, q)

			if i >= p && i < q {
				i -= p
				fmt.Printf("zeroing %d %d\n", i, p)
				y[yPtr] = 0
				yPtr++
			} else {
				s := int32(0)
				if i >= q {
					i -= q
					s = -1
				}

				k0 := k
				for {
					k--
					rowK := pvqURow(int(k))
					pVal := uint32(rowK[n])
					if i >= pVal {
						break
					}
				}
				i -= p
				fmt.Printf("i %d p %d\n", i, p)
				y[yPtr] = update(k0, k, s, &norm)
				yPtr++
			}
		}
		n--
	}

	// 处理 n=2
	{
		p := uint32(2*k + 1)
		fmt.Printf("p %d i %d\n", p, i)
		s := int32(0)
		if i >= p {
			i -= p
			s = -1
		}

		k0 := k
		k = (i + 1) / 2
		if k != 0 {
			i -= 2*k - 1
		}

		fmt.Println("n == 2")
		y[yPtr] = update(k0, k, s, &norm)
		yPtr++
	}

	// 处理 n=1
	{
		s := int32(-int64(i))
		fmt.Println("n == 1")
		y[yPtr] = update(k, 0, s, &norm)
	}

	fmt.Printf("norm %d\n", norm)
	return norm
}

func pvqU(n, k int) int {
	min := n
	if k < n {
		min = k
	}
	max := n
	if k > n {
		max = k
	}
	return int(pvqURow(min)[max])
}

func pvqV(n, k int) int {
	return pvqU(n, k) + pvqU(n, k+1)
}

func decodePulses(rd *comm.RangeDecoder, y []int32, n, k int) float32 {
	idx := rd.DecodeUniform(pvqV(n, k))
	fmt.Printf("idx %d\n", idx)
	return float32(cwrsi(uint32(n), uint32(k), uint32(idx), y))
}

// expRotation1 实现点旋转操作
func expRotation1(x []float32, len, stride int, c, s float32) {
	end := len - stride
	for i := 0; i < end; i++ {
		x1 := x[i]
		x2 := x[i+stride]

		x[i+stride] = c*x2 + s*x1
		x[i] = c*x1 - s*x2
	}

	// 反向处理
	for i := end - stride - 1; i >= 0; i-- {
		x1 := x[i]
		x2 := x[i+stride]
		x[i+stride] = c*x2 + s*x1
		x[0] = c*x1 - s*x2 // 固定更新第一个元素
	}
}

// ExpRotation 实现指数旋转
func ExpRotation(x []float32, _len, stride, k, spread int) {
	if 2*k >= _len || spread == SPREAD_NONE {
		return
	}

	gain := float32(_len) / float32(_len+(20-5*spread)*k)
	theta := float32(PI) * gain * gain / 4.0

	c := float32(math.Cos(float64(theta)))
	s := float32(math.Sin(float64(theta)))

	stride2 := 0
	if _len >= stride<<3 {
		stride2 = 1
		// 计算等效于 sqrt(len/stride) 的整数值
		for (stride2*stride2+stride2)*stride+(stride>>2) < _len {
			stride2++
		}
	}

	for i := 0; i < stride; i++ {
		start := i * _len
		if start+_len > len(x) { // 防止越界
			continue
		}
		if stride2 != 0 {
			expRotation1(x[start:], _len, stride2, s, c)
		}
		expRotation1(x[start:], _len, 1, c, s)
	}
}

// ExtractCollapseMask 提取塌陷掩码
func ExtractCollapseMask(y []int32, b int) int {
	if b <= 1 {
		return 1
	}

	collapseMask := 0
	for i := 0; i < len(y); i += b {
		for j := 0; j < b; j++ {
			if i+j >= len(y) {
				break
			}
			if y[i+j] != 0 {
				collapseMask |= 1 << j
			}
		}
	}

	return collapseMask
}

// Bits2Pulses 将比特数转换为脉冲数
func Bits2Pulses(cache []uint8, bits int32) int32 {
	low := 0
	high := len(cache) - 1

	for i := 0; i < 6; i++ { // 二分查找
		center := (low + high + 1) >> 1
		if int32(cache[center]) >= bits {
			high = center
		} else {
			low = center
		}
	}

	lowBits := bits - int32(cache[low])
	if low == 0 {
		lowBits = bits + 1 // 特殊处理第一个元素
	}
	highBits := int32(cache[high]) - bits

	if lowBits <= highBits {
		return int32(low)
	}
	return int32(high)
}

// Pulses2Bits 将脉冲数转换为比特数
func Pulses2Bits(cache []uint8, pulses int32) int32 {
	if pulses == 0 || int(pulses) >= len(cache) {
		return 0
	}
	return int32(cache[pulses])
}

// Unquantize 实现去量化
func Unquantize(
	rd *comm.RangeDecoder,
	x []float32,
	n, k, spread, blocks int,
	gain float32,
) int {
	y := make([]int32, 176) // 显式初始化，避免 unsafe

	// 解码脉冲并计算增益
	pulseNorm := float32(decodePulses(rd, y, n, k))
	gain /= float32(math.Sqrt(float64(pulseNorm)))

	// 应用增益
	for i := 0; i < n; i++ {
		x[i] = gain * float32(y[i])
	}

	// 应用指数旋转
	ExpRotation(x, n, blocks, k, spread)

	// 提取塌陷掩码
	return ExtractCollapseMask(y, blocks)
}

// RenormalizeVector 向量重新归一化
func RenormalizeVector(x []float32, gain float32) {
	g := float32(0.0)
	for _, v := range x {
		g += v * v
	}

	if g == 0 {
		return // 防止除以零
	}
	gain /= float32(math.Sqrt(float64(g)))

	for i := range x {
		x[i] *= gain
	}
}

// StereoMerge 立体声合并
func StereoMerge(x, y []float32, mid float32, n int) {
	xp := float32(0.0)
	side := float32(0.0)

	// 计算点积和边信息
	for i := 0; i < n; i++ {
		xp += x[i] * y[i]
		side += y[i] * y[i]
	}

	xp *= mid // 调整中间值

	e := mid*mid + side

	e0 := e - 2*xp
	e1 := e + 2*xp

	if e0 < 6e-4 || e1 < 6e-4 {
		copy(y[:n], x[:n]) // 直接复制
		return
	}

	gain0 := float32(1.0) / float32(math.Sqrt(float64(e0)))
	gain1 := float32(1.0) / float32(math.Sqrt(float64(e1)))

	// 应用立体声合并
	for i := 0; i < n; i++ {
		v0 := mid * x[i]
		v1 := y[i]
		x[i] = gain0 * (v0 - v1)
		y[i] = gain1 * (v0 + v1)
	}
}

type BandInfo struct {
	Itheta int     // 原 itheta: usize
	Inv    bool    // 原 inv: bool
	Mid    float32 // 原 mid: f32
	Side   float32 // 原 side: f32
	Delta  int16   // 原 delta: i16
	Qalloc int     // 原 qalloc: usize
	Fill   int     // 原 fill: usize
}

// NewCelt 构造函数
func NewCelt(stereo bool) *Celt {
	return &Celt{
		stereo:           stereo,
		stereo_pkt:       false,
		bits:             0,
		lm:               0,
		band_start:       0,
		band_end:         MAX_BANDS,
		spread:           SPREAD_NORMAL,
		anticollapse_bit: 0,
		blocks:           0,
		blocksize:        0,
	}
}

// Setup 方法
func (c *Celt) Setup(pkt struct{ stereo bool }) {
	c.stereo_pkt = pkt.stereo
}

// ResetGains 方法
func (c *Celt) ResetGains() {
	for i := range c.frames {
		c.frames[i].pf.gains_new = [3]float32{0, 0, 0}
	}
}

// ParsePostfilter 方法
func (c *Celt) ParsePostfilter(rd *comm.RangeDecoder) {
	if rd.DecodeLogP(1) {
		octave := rd.DecodeUniform(6)
		period := (16 << octave) + int(rd.RawBits(4+octave)) - 1
		gainBits := int(rd.RawBits(3)) + 1
		gain := float32(gainBits) * 0.09375

		tapset := 0
		if rd.Available() >= 2 {
			tapset = rd.DecodeICDF(TAPSET)
		}

		fmt.Printf(
			"postfilter: octave %d, period %d, gain %f, tapset %d\n",
			octave, period, gain, tapset,
		)

		taps := POSTFILTER_TAPS[tapset]
		for i := range c.frames {
			if period < MIN_PERIOD {
				c.frames[i].pf.period_new = MIN_PERIOD
			} else {
				c.frames[i].pf.period_new = period
			}
			c.frames[i].pf.gains_new[0] = taps[0] * gain
			c.frames[i].pf.gains_new[1] = taps[1] * gain
			c.frames[i].pf.gains_new[2] = taps[2] * gain
		}
	} else {
		fmt.Println("postfilter: no")
	}
}

// DecodeCoarseEnergy 方法
func (c *Celt) DecodeCoarseEnergy(rd *comm.RangeDecoder, bandStart, bandEnd int) {
	var alpha, beta float32
	var model []uint8

	if rd.Available() > 3 && rd.DecodeLogP(3) {
		alpha = 0.0
		beta = 1.0 - 4915.0/32768.0
		model = COARSE_ENERGY_INTRA[c.lm]
	} else {
		alpha = ALPHA_COEF[c.lm]
		beta = BETA_COEF[c.lm]
		model = COARSE_ENERGY_INTER[c.lm]

	}

	fmt.Printf("model %.6f %.6f\n", alpha, beta)

	prev := [MAX_FRAMES]float32{0, 0}

	for i := 0; i < MAX_BANDS; i++ {
		for j := 0; j < len(c.frames); j++ {
			frame := &c.frames[j]
			en := &frame.energy[i]

			if i < bandStart || i >= bandEnd {
				*en = 0.0
				continue
			}

			available := rd.Available()
			fmt.Printf("available %d\n", available)

			var value float32
			if available >= 15 {
				k := i
				if k > 20 {
					k = 20
				}
				k <<= 1 // k * 2

				v := float32(rd.DecodeLaplace(int(model[k]<<7), int(model[k+1]<<6)))
				fmt.Printf("decode_laplace %.6f <- %d %d\n", v, i, k)
				value = v
			} else if available >= 1 {
				v := rd.DecodeICDF(MODEL_ENERGY_SMALL)
				// (v >> 1) ^ -(v & 1) 的等效计算
				value = float32(v >> 1)
				if (v & 1) != 0 {
					value = -value - 1
				}
			} else {
				value = -1
			}

			fmt.Printf(
				"energy %d/%d %.6f * %.6f + %.6f + %.6f\n",
				i, j, *en, alpha, prev[j], value,
			)

			// 更新能量值 *en = max(en, -9) * alpha + prev[j] + value
			current := *en
			if current < -9.0 {
				current = -9.0
			}
			newEn := current*alpha + prev[j] + value

			// 防止NaN或无限值
			if math.IsNaN(float64(newEn)) || math.IsInf(float64(newEn), 0) {
				newEn = 0.0
			}

			*en = newEn
			prev[j] += beta * value

			// 限制prev[j]的范围防止溢出
			if prev[j] > 1000 {
				prev[j] = 1000
			} else if prev[j] < -1000 {
				prev[j] = -1000
			}
		}

		// 如果是单声道包，只处理第一个frame
		if !c.stereo_pkt {
			break
		}
	}

	fmt.Printf("Frame0 energy: %.6v\n", c.frames[0].energy[:])
	if c.stereo_pkt {
		fmt.Printf("Frame1 energy: %.6v\n", c.frames[1].energy[:])
	}
}

// DecodeTfChanges 方法
func (c *Celt) DecodeTfChanges(rd *comm.RangeDecoder, bandStart, bandEnd int, transient bool) {
	tfChanged := make([]bool, MAX_BANDS)

	bits := struct{ field0, field1 int }{2, 4}
	if !transient {
		bits = struct{ field0, field1 int }{4, 5}
	}

	available := rd.Available()
	tfSelect := TF_SELECT[c.lm][boolToInt(transient)]

	selectBit := (c.lm != 0) && (available > bits.field0)
	fmt.Printf("selectBit %t %d\n", selectBit, available)

	fieldBits := bits.field0
	diff := false
	changed := false

	for i := bandStart; i < bandEnd; i++ {
		if available > fieldBits+boolToInt(selectBit) {
			// 模拟 diff ^= rd.decode_logp(field_bits)
			if rd.DecodeLogP(fieldBits) {
				diff = !diff
			}
			fmt.Printf("band %d bits %d %t\n", i, fieldBits, diff)
			available = rd.Available()
			changed = changed || diff
		}

		tfChanged[i] = diff
		fieldBits = bits.field1
	}

	selectIdx := 0
	if selectBit {
		changedInt := boolToInt(changed)
		if tfSelect[0][changedInt] != tfSelect[1][changedInt] {
			if rd.DecodeLogP(1) {
				selectIdx = 1
			}
		}
	}

	for i := bandStart; i < bandEnd; i++ {
		c.tf_change[i] = tfSelect[selectIdx][boolToInt(tfChanged[i])]
	}

	fmt.Printf("tf_change %#v\n", c.tf_change[bandStart:bandEnd])
}

// 辅助函数：bool转int
func boolToInt(b bool) int {
	if b {
		return 1
	}
	return 0
}

func (d *Celt) DecodeAllocation(rd *comm.RangeDecoder, bandStart, bandEnd int) {
	band := struct{ start, end int }{bandStart, bandEnd}
	var caps [MAX_BANDS]int32
	var threshold [MAX_BANDS]int32
	var trim_offset [MAX_BANDS]int32
	var boost [MAX_BANDS]int32

	scale := d.lm
	if d.stereo_pkt {
		scale++
	}

	spread := SPREAD_NORMAL
	if rd.Available() > 4 {
		spread = rd.DecodeICDF(MODEL_SPREAD)
	}
	_ = spread // 根据实际需要使用

	staticCaps := STATIC_CAPS[d.lm][boolToInt(d.stereo_pkt)]

	for i := 0; i < MAX_BANDS; i++ {
		if i < len(staticCaps) && i < len(FREQ_RANGE) {
			caps[i] = int32((staticCaps[i] + 64) * FREQ_RANGE[i] << scale >> 2)
		}
	}

	dynalloc := 6
	boostSize := 0

	for i := band.start; i < band.end; i++ {
		quanta := FREQ_RANGE[i] << scale
		if quanta<<3 > quanta {
			if quanta > 6<<3 {
				quanta = 6 << 3
			} else {
				quanta = quanta << 3
			}
		}
		bandDynalloc := dynalloc
		for (bandDynalloc<<3)+boostSize < rd.AvailableFrac() && boost[i] < caps[i] {
			add := rd.DecodeLogP(bandDynalloc)
			if !add {
				break
			}
			boost[i] += int32(quanta)
			boostSize += int(quanta)
			bandDynalloc = 1
		}

		if boost[i] != 0 && dynalloc > 2 {
			dynalloc--
		}
	}

	allocTrim := 5
	if rd.AvailableFrac() > boostSize+(6<<3) {
		allocTrim = rd.DecodeICDF(ALLOC_TRIM)
	}

	available := rd.AvailableFrac() - 1
	d.anticollapse_bit = 0
	if d.blocks > 1 && d.lm >= 2 && available >= (d.lm+2)<<3 {
		available -= 1 << 3
		d.anticollapse_bit = 1 << 3
	}

	skipBit := 0
	if available >= 1<<3 {
		available -= 1 << 3
		skipBit = 1 << 3
	}

	intensityStereoBit := 0
	dualStereoBit := 0
	if d.stereo_pkt {
		intensityStereo := int(LOG2_FRAC[band.end-band.start])
		if intensityStereo <= available {
			available -= intensityStereo
			intensityStereoBit = intensityStereo
			if available >= 1<<3 {
				available -= 1 << 3
				dualStereoBit = 1 << 3
			}
		}
	}

	for i := band.start; i < band.end; i++ {
		trim := allocTrim - (5 + d.lm)
		rangeVal := int32(FREQ_RANGE[i]) * int32(band.end-i-1)
		lm := d.lm + 3
		stereoScale := lm
		if d.stereo_pkt {
			stereoScale++
		}
		var stereoThreshold int32 = 0
		if d.stereo_pkt {
			stereoThreshold = 1 << 8
		}

		threshold[i] = int32((3 * FREQ_RANGE[i] << lm) >> 4)
		if threshold[i] < stereoThreshold {
			threshold[i] = stereoThreshold
		}

		trim_offset[i] = int32(trim) * (rangeVal << stereoScale) >> 6

		if FREQ_RANGE[i]<<d.lm == 1 {
			trim_offset[i] -= stereoThreshold
		}
	}

	var codedChannelBits int32 = 1
	if d.stereo_pkt {
		codedChannelBits = 2
	}
	codedChannelBits <<= 3

	low := 1
	high := CELT_VECTOR - 1
	for low <= high {
		center := (low + high) / 2
		done := false
		var total int32 = 0

		for i := band.end - 1; i >= band.start; i-- {
			bandBits := int32(FREQ_RANGE[i] * STATIC_ALLOC[center][i])
			if d.stereo_pkt {
				bandBits <<= 1
			}
			bandBits <<= d.lm
			bandBits >>= 2

			if bandBits != 0 {
				bandBits += trim_offset[i]
				if bandBits < 0 {
					bandBits = 0
				}
			}

			bandBits += boost[i]

			if bandBits >= threshold[i] || done {
				done = true
				if bandBits > caps[i] {
					total += caps[i]
				} else {
					total += bandBits
				}
			} else {
				if bandBits >= codedChannelBits {
					total += codedChannelBits
				}
			}
		}

		if total > int32(available) {
			high = center - 1
		} else {
			low = center + 1
		}
	}

	high = low
	low = high - 1

	var bits1 [MAX_BANDS]int32
	var bits2 [MAX_BANDS]int32
	skipStartband := band.start

	bitsEstimation := func(idx, i int) int32 {
		bits := int32(FREQ_RANGE[i] * STATIC_ALLOC[idx][i])
		if d.stereo_pkt {
			bits <<= 1
		}
		bits <<= d.lm
		bits >>= 2

		if bits != 0 {
			bits += trim_offset[i]
			if bits < 0 {
				bits = 0
			}
		}
		return bits
	}

	for i := band.start; i < band.end; i++ {
		bits1[i] = bitsEstimation(low, i)
		bits2[i] = bitsEstimation(high, i)

		if boost[i] != 0 {
			if low != 0 {
				bits1[i] += boost[i]
			}
			bits2[i] += boost[i]
			skipStartband = i
		}

		bits2[i] -= bits1[i]
		if bits2[i] < 0 {
			bits2[i] = 0
		}
	}

	low = 0
	high = 1 << ALLOC_STEPS

	for step := 0; step < ALLOC_STEPS; step++ {
		center := int32((low + high) / 2)
		done := false
		var total int32 = 0

		for j := band.end - 1; j >= band.start; j-- {
			bits := bits1[j] + (center * bits2[j] >> ALLOC_STEPS)

			if bits >= threshold[j] || done {
				done = true
				if bits > caps[j] {
					total += caps[j]
				} else {
					total += bits
				}
			} else if bits >= codedChannelBits {
				total += codedChannelBits
			}
		}

		if total > int32(available) {
			high = int(center)
		} else {
			low = int(center)
		}
	}

	done := false
	total := 0
	for i := band.end - 1; i >= band.start; i-- {
		bits := bits1[i] + (int32(low) * bits2[i] >> ALLOC_STEPS)

		if bits < threshold[i] && !done {
			if bits >= codedChannelBits {
				bits = codedChannelBits
			} else {
				bits = 0
			}
		} else {
			done = true
		}

		if bits > caps[i] {
			bits = caps[i]
		}
		d.pulses[i] = bits
		total += int(bits)
	}

	codedband := band.end
bandsLoop:
	for j := band.end - 1; j >= band.start; j-- {
		if j == skipStartband {
			available += skipBit
			codedband = j + 1
			break bandsLoop
		}

		bandDelta := int32(FREQ_BANDS[j+1] - FREQ_BANDS[band.start])
		remaining := available - total
		bits := int32(remaining) / bandDelta
		remaining -= int(bits * bandDelta)

		allocation := d.pulses[j] + bits*int32(FREQ_BANDS[j])
		if remaining-int(bandDelta) > 0 {
			allocation += int32(remaining) - bandDelta
		}

		minVal := codedChannelBits
		if threshold[j] > minVal {
			minVal = threshold[j]
		}
		if allocation >= minVal {
			if rd.DecodeLogP(1) {
				codedband = j + 1
				break bandsLoop
			}
			total += 1 << 3
			allocation -= 1 << 3
		}

		total -= d.pulses[j]
		if intensityStereoBit != 0 {
			total -= intensityStereoBit
			intensityStereoBit = LOG2_FRAC[j-band.start]
			total += intensityStereoBit
		}

		if allocation >= codedChannelBits {
			d.pulses[j] = codedChannelBits
		} else {
			d.pulses[j] = 0
		}
		total += d.pulses[j]
	}

	d.intensity_stereo = 0
	if intensityStereoBit != 0 {
		d.intensity_stereo = band.start + rd.DecodeUniform(codedband-band.start)
	}

	d.dual_stereo = false
	if d.intensity_stereo <= band.start {
		available += dualStereoBit
	} else if dualStereoBit != 0 {
		d.dual_stereo = rd.DecodeLogP(1)
	}

	bandDelta := FREQ_BANDS[codedband] - FREQ_BANDS[band.start]
	remaining := available - total
	bandbits := remaining / bandDelta
	remaining -= bandbits * bandDelta

	for i := band.start; i < band.end; i++ {
		fr := FREQ_RANGE[i]
		bits := bandbits * fr
		if remaining > 0 {
			if remaining < fr {
				bits += remaining
				remaining = 0
			} else {
				bits += fr
				remaining -= fr
			}
		}
		d.pulses[i] += bits
	}

	extrabits := 0
	for i := band.start; i < band.end; i++ {
		n := FREQ_RANGE[i] << d.lm
		prevExtra := extrabits
		d.pulses[i] += extrabits

		if n > 1 {
			extrabits = d.pulses[i] - caps[i]
			if extrabits < 0 {
				extrabits = 0
			}
			d.pulses[i] -= extrabits

			dof := n
			if d.stereo_pkt {
				dof *= 2
				if n > 2 && !d.dual_stereo && i < d.intensity_stereo {
					dof--
				}
			}

			duration := d.lm << 3
			dofChannels := dof * (LOG_FREQ_RANGE[i] + duration)
			offset := (dofChannels >> 1) - dof*FINE_OFFSET

			if n == 2 {
				offset += dof * 2
			}

			pulse := d.pulses[i] + offset
			if pulse < 2*(dof<<3) {
				offset += dofChannels >> 2
			} else if pulse < 3*(dof<<3) {
				offset += dofChannels >> 3
			}

			pulse = d.pulses[i] + offset
			fineBits := (pulse + (dof << 2)) / (dof << 3)
			maxBits := d.pulses[i] >> 3
			if d.stereo_pkt {
				maxBits >>= 1
			}
			if maxBits > MAX_FINE_BITS {
				maxBits = MAX_FINE_BITS
			}
			if fineBits < 0 {
				fineBits = 0
			} else if fineBits > maxBits {
				fineBits = maxBits
			}
			d.fine_bits[i] = fineBits
			d.fine_priority[i] = fineBits*(dof<<3) >= pulse

			d.pulses[i] -= fineBits << boolToInt(d.stereo_pkt) << 3
		} else {
			extrabits = d.pulses[i] - boolToInt(d.stereo_pkt)*8 - 8
			if extrabits < 0 {
				extrabits = 0
			}
			d.pulses[i] -= extrabits
			d.fine_bits[i] = 0
			d.fine_priority[i] = true
		}

		if extrabits > 0 {
			scale := boolToInt(d.stereo_pkt) + 3
			extraFine := MAX_FINE_BITS - d.fine_bits[i]
			if extraFine > extrabits>>scale {
				extraFine = extrabits >> scale
			}
			d.fine_bits[i] += extraFine
			extraFineBits := extraFine << scale
			d.fine_priority[i] = extraFineBits >= extrabits-prevExtra
			extrabits -= extraFineBits
		}
	}

	d.remaining = extrabits

	for i := codedband; i < band.end; i++ {
		d.fine_bits[i] = d.pulses[i] >> boolToInt(d.stereo_pkt) >> 3
		d.pulses[i] = 0
		d.fine_priority[i] = d.fine_bits[i] < 1
	}

	d.codedband = codedband
}

func (c *Celt) DecodeFineEnergy(rd *comm.RangeDecoder, bandStart, bandEnd int) {
	for i := bandStart; i < bandEnd; i++ {
		fineBits := c.fine_bits[i]
		if fineBits == 0 {
			continue
		}

		for f := 0; f <= boolToInt(c.stereo_pkt); f++ {
			frame := &c.frames[f]
			q2 := float32(rd.RawBits(fineBits))
			fmt.Printf("-- fine_bits %d\n", fineBits)

			offset := (q2+0.5)*float32(1<<(14-fineBits))/16384.0 - 0.5
			fmt.Printf("q2 %.6f offset %.6f\n", q2, offset)

			frame.energy[i] += offset
		}
	}
}

func (c *Celt) DecodeBand1(rd *comm.RangeDecoder, midBuf []float32, sideBuf, lowbandOut []float32) {
	oneSample := func() float32 {
		if c.remaining2 >= 1<<3 {
			c.remaining2 -= 1 << 3
			if rd.RawBits(1) != 0 {
				return -1.0
			}
		}
		return 1.0
	}

	midBuf[0] = oneSample()
	if sideBuf != nil {
		sideBuf[0] = oneSample()
	}
	if lowbandOut != nil {
		lowbandOut[0] = midBuf[0]
	}
}

func (c *Celt) ComputeTheta(
	rd *comm.RangeDecoder,
	band, lm, n, b, b0, blocks int,
	dualstereo bool,
	fill int,
) BandInfo {
	fmt.Printf("band %d\n", band)
	pulseCap := LOG_FREQ_RANGE[band] + lm*8
	offset := pulseCap / 2
	if dualstereo && n == 2 {
		offset -= QTHETA_OFFSET_TWOPHASE
	} else {
		offset -= QTHETA_OFFSET
	}

	qn := 1
	if !(dualstereo && band >= c.intensity_stereo) {
		n2 := 2*n - 1
		if dualstereo && n == 2 {
			n2 = 2*n - 2
		}
		fmt.Printf("n2 %d pulse_cap %d b %d\n", n2, pulseCap, b)

		qb := b - pulseCap - (4 << 3)
		if qb > (b+n2*offset)/n2 {
			qb = (b + n2*offset) / n2
		}
		if qb > 8<<3 {
			qb = 8 << 3
		}

		if qb >= (1<<3)/2 {
			idx := qb & 0x7
			shift := 14 - (qb >> 3)
			qn = (int(QN_EXP2[idx]) >> shift) + 1
			qn = (qn >> 1) << 1
		}
	}

	fmt.Printf("qn %d\n", qn)

	tellFrac := rd.TellFrac()
	itheta := 0
	inv := false

	if qn != 1 {
		if dualstereo && n > 2 {
			itheta = rd.DecodeStep(qn / 2)
		} else if dualstereo || b0 > 1 {
			itheta = rd.DecodeUniform(qn + 1)
		} else {
			itheta = rd.DecodeTriangular(qn)
		}
		itheta = itheta * 16384 / qn
	} else if dualstereo && b > BITRES && c.remaining2 > BITRES {
		inv = rd.DecodeLogP(2)
	}

	qalloc := rd.TellFrac() - tellFrac

	imid := 0
	iside := 0
	newFill := fill
	delta := 0

	switch itheta {
	case 0:
		imid = 32767
		iside = 0
		newFill = fill & ((1 << blocks) - 1)
		delta = -16384
	case 16384:
		imid = 0
		iside = 32767
		newFill = fill & (((1 << blocks) - 1) << blocks)
		delta = 16384
	default:
		imid = int(bitexactCos(int16(itheta)))
		iside = int(bitexactCos(int16(16384 - itheta)))
		log2tan := bitexactLog2tan(int16(iside), int16(imid))
		delta = int(bitexactFracMul16(int16((n-1)<<7), int16(log2tan)))
	}

	return BandInfo{
		Itheta: itheta,
		Inv:    inv,
		Mid:    float32(imid) / 32768.0,
		Side:   float32(iside) / 32768.0,
		Delta:  int16(delta),
		Qalloc: qalloc,
		Fill:   newFill,
	}
}

func (c *Celt) Rng() uint32 {
	c.seed = 1664525*c.seed + 1013904223
	return c.seed
}

func (c *Celt) DecodeBandNoSplit(
	rd *comm.RangeDecoder,
	midBuf []float32,
	lowband []float32,
	n, blocks int,
	gain float32,
	cache []uint8,
	b int,
	fill int,
) int {
	q := Bits2Pulses(cache, int32(b))
	currBits := Pulses2Bits(cache, q)

	c.remaining2 -= int32(currBits)

	for c.remaining2 < 0 && q > 0 {
		c.remaining2 += int32(currBits)
		q--
		currBits = Pulses2Bits(cache, q)
		c.remaining2 -= int32(currBits)
	}

	if q != 0 {
		k := q
		if q >= 8 {
			k = (8 + (q & 7)) << ((q >> 3) - 1)
		}
		return Unquantize(rd, midBuf, n, int(k), c.spread, blocks, gain)
	}

	cmMask := (1 << blocks) - 1
	newFill := fill & cmMask
	if newFill == 0 {
		for i := 0; i < n; i++ {
			midBuf[i] = 0
		}
		return 0
	}

	cm := cmMask
	if lowband != nil {
		for i := 0; i < n; i++ {
			rng := float32(0)
			if (c.Rng() & 0x8000) != 0 {
				rng = 1.0 / 256.0
			} else {
				rng = -1.0 / 256.0
			}
			midBuf[i] = lowband[i] + rng
		}
	} else {
		for i := 0; i < n; i++ {
			midBuf[i] = float32(int32(c.Rng()) >> 20)
		}
	}

	RenormalizeVector(midBuf[:n], gain)
	return cm
}

func (c *Celt) DecodeBand(
	rd *comm.RangeDecoder,
	band int,
	midBuf []float32,
	sideBuf []float32,
	n int,
	b int,
	blocks int,
	lowband []float32,
	lowbandOut []float32,
	lm int,
	level int,
	gain float32,
	fill int,
) int {
	nB := n / blocks
	nB0 := nB
	dualstereo := sideBuf != nil
	n0 := n
	b0 := blocks

	timeDivide := 0
	longblocks := b0 == 1
	fmt.Printf(
		"decode_band N=%d lowband_out %t\n",
		n, lowbandOut != nil,
	)

	fmt.Println("mid_buf")
	for i, v := range midBuf[:n] {
		fmt.Printf("%d: %.08f\n", i, v)
	}

	if n == 1 {
		c.DecodeBand1(rd, midBuf, sideBuf, lowbandOut)
		return 1
	}

	var lowbandCopy []float32
	if lowband != nil {
		lowbandCopy = make([]float32, n)
		copy(lowbandCopy, lowband[:n])
	}

	recombine := 0
	if !dualstereo && level == 0 {
		tfChange := c.tf_change[band]
		if tfChange < 0 {
			recombine = 0
		} else {
			recombine = tfChange
		}

		fmt.Printf("recombine %d\n", recombine)

		if lowbandCopy != nil {
			fmt.Println("lowband")
			for i, v := range lowbandCopy {
				fmt.Printf("%d: %.08f\n", i, v)
			}
		}

		for k := 0; k < recombine; k++ {
			if lowbandCopy != nil {
				haar1(lowbandCopy, n>>k, 1<<k)
			}

			fill = int(BIT_INTERLEAVE[fill&0xf] | (BIT_INTERLEAVE[fill>>4] << 2))
		}

		blocks >>= recombine
		nB <<= recombine
		fmt.Printf("blocks %d N_B %d\n", blocks, nB)
		for (nB&1) == 0 && tfChange < 0 {
			if lowbandCopy != nil {
				fmt.Println("EDIT")
				haar1(lowbandCopy, nB, blocks)
			}

			fill |= fill << blocks
			blocks <<= 1
			nB >>= 1

			timeDivide++
			tfChange++
		}

		b0 = blocks
		nB0 = nB

		fmt.Printf("B0 %d\n", b0)
		if b0 > 1 && lowbandCopy != nil {
			deinterleaveHadamard(
				c.scratch[:],
				lowbandCopy,
				nB>>recombine,
				b0<<recombine,
				longblocks,
			)
		}
	}

	lowband = lowbandCopy

	cacheStart := CACHE_INDEX[(lm+1)*MAX_BANDS+band]
	cache := CACHE_BITS[cacheStart:]

	split := false
	if !dualstereo && lm >= 0 && b > (int(cache[0])+12) && n > 2 {
		n >>= 1
		splitMid := midBuf[:n]
		splitSide := midBuf[n : n*2]
		midBuf = splitMid
		sideBuf = splitSide
		lm--
		if blocks == 1 {
			fill = (fill & 1) | (fill << 1)
		}
		blocks = (blocks + 1) >> 1
		split = true
	} else {
		split = dualstereo
	}

	fmt.Printf("split %t blocks %d lm %d\n", split, blocks, lm)

	cm := 0
	if sideBuf != nil {
		info := c.ComputeTheta(
			rd,
			band,
			lm,
			n,
			b,
			b0,
			blocks,
			dualstereo,
			fill,
		)
		itheta := info.Itheta
		inv := info.Inv
		mid := info.Mid
		side := info.Side
		delta := info.Delta
		qalloc := info.Qalloc

		b -= qalloc

		fmt.Printf(
			"itheta %d delta %d n %d dualstereo %t\n",
			itheta, delta, n, dualstereo,
		)

		if n == 2 && dualstereo {
			sbits := 0
			if itheta != 0 && itheta != 16384 {
				sbits = 1 << 3
			}

			mbits := b - sbits
			c.remaining2 -= qalloc + sbits

			var midBuf2, sideBuf2 []float32
			if itheta > 8192 {
				midBuf2 = sideBuf[:n]
				sideBuf2 = midBuf[:n]
			} else {
				midBuf2 = midBuf[:n]
				sideBuf2 = sideBuf[:n]
			}

			sign := float32(1.0)
			if sbits != 0 {
				if rd.RawBits(1) != 0 {
					sign = -1.0
				}
			}

			cm = c.DecodeBand(
				rd,
				band,
				midBuf2,
				nil,
				n,
				mbits,
				blocks,
				lowband,
				lowbandOut,
				lm,
				level,
				gain,
				fill,
			)
			sideBuf2[1] = -sign * midBuf2[0]
		} else {
			nextLowband2 := []float32(nil)
			nextLowbandOut1 := []float32(nil)

			if b0 > 1 && !dualstereo && (itheta&0x3fff) != 0 {
				if itheta > 8192 {
					delta -= delta >> (4 - lm)
				} else {
					newDelta := delta + (n<<3)>>(5-lm)
					if newDelta < 0 {
						delta = newDelta
					}
				}
			}

			fmt.Printf("delta %d\n", delta)
			mbits := (b - delta) / 2
			if mbits < 0 {
				mbits = 0
			} else if mbits > b {
				mbits = b
			}
			sbits := b - mbits

			c.remaining2 -= qalloc

			if !dualstereo && lowband != nil {
				nextLowband2 = lowband[n:]
				lowband = lowband[:n]
			}

			nextLevel := level
			if !dualstereo {
				nextLevel++
			}

			rebalance := c.remaining2
			sideShift := 0
			if dualstereo {
				sideShift = b0 >> 1
			}

			if mbits >= sbits {
				cm = c.DecodeBand(
					rd,
					band,
					midBuf,
					nil,
					n,
					mbits,
					blocks,
					lowband,
					nextLowbandOut1,
					lm,
					nextLevel,
					gain*mid,
					fill,
				)

				rebalance = mbits - (rebalance - c.remaining2)
				if rebalance > 3<<3 && itheta != 0 {
					sbits += rebalance - (3 << 3)
				}

				sideCM := c.DecodeBand(
					rd,
					band,
					sideBuf,
					nil,
					n,
					sbits,
					blocks,
					nextLowband2,
					nil,
					lm,
					nextLevel,
					gain*side,
					fill>>blocks,
				)
				cm |= sideCM << sideShift
			} else {
				sideCM := c.DecodeBand(
					rd,
					band,
					sideBuf,
					nil,
					n,
					sbits,
					blocks,
					nextLowband2,
					nil,
					lm,
					nextLevel,
					gain*side,
					fill>>blocks,
				)
				cm = sideCM << sideShift

				rebalance = sbits - (rebalance - c.remaining2)
				if rebalance > 3<<3 && itheta != 16384 {
					mbits += rebalance - (3 << 3)
				}

				midCM := c.DecodeBand(
					rd,
					band,
					midBuf,
					nil,
					n,
					mbits,
					blocks,
					lowband,
					nextLowbandOut1,
					lm,
					nextLevel,
					gain*mid,
					fill,
				)
				cm |= midCM
			}
		}
	} else {
		cm = c.DecodeBandNoSplit(
			rd,
			midBuf,
			lowband,
			n,
			blocks,
			gain,
			cache,
			b,
			fill,
		)
	}

	if level == 0 {
		if b0 > 1 {
			InterleaveHadamard(
				c.scratch[:],
				midBuf,
				nB>>recombine,
				b0<<recombine,
				longblocks,
			)
		}

		nB = nB0
		blocks = b0
		for i := 0; i < timeDivide; i++ {
			blocks >>= 1
			nB <<= 1
			cm |= cm >> blocks
			Haar1(midBuf, nB, blocks)
		}

		for k := 0; k < recombine; k++ {
			cm = BIT_DEINTERLEAVE[cm]
			Haar1(midBuf, n0>>k, 1<<k)
		}

		blocks <<= recombine

		if lowbandOut != nil {
			nSqrt := float32(math.Sqrt(float64(n0)))
			fmt.Println("Lowband_out")
			for i := 0; i < n0; i++ {
				lowbandOut[i] = nSqrt * midBuf[i]
				fmt.Printf("%.08f\n", lowbandOut[i])
			}
		}

		cm &= (1 << blocks) - 1
		fmt.Printf("cm %d\n", cm)
	}

	return cm
}

func (c *Celt) DecodeBands(
	rd *comm.RangeDecoder,
	bandStart, bandEnd int,
	coeff0, coeff1 []float32,
) {
	lm := c.lm
	updateLowband := true
	lowbandOffset := 0

	const NORM_SIZE = 8 * 100
	normMid := make([]float32, NORM_SIZE)
	normSide := make([]float32, NORM_SIZE)

	for i := bandStart; i < bandEnd; i++ {
		bandOffset := int(FREQ_BANDS[i] << lm)
		bandSize := int(FREQ_RANGE[i] << lm)

		x := coeff0[bandOffset : bandOffset+bandSize]
		y := coeff1[bandOffset : bandOffset+bandSize]

		consumed := rd.TellFrac()

		if i != bandStart {
			c.remaining -= consumed
		}

		c.remaining2 = rd.AvailableFrac() - 1 - c.anticollapse_bit

		b := 0
		if i <= c.codedband-1 {
			fmt.Printf("rem %d rem2 %d\n", c.remaining, c.remaining2)
			remaining := c.remaining / min(c.codedband-1, 3)
			b = min(c.remaining2+1, c.pulses[i]+remaining)
			b = max(0, b)
			b = min(16383, b)
		}

		fmt.Printf("b %d\n", b)

		if FREQ_BANDS[i]-FREQ_RANGE[i] >= FREQ_BANDS[bandStart] &&
			(updateLowband || lowbandOffset == 0) {
			lowbandOffset = i
		}

		cm := [2]int{0, 0}
		var effectiveLowband int
		if lowbandOffset != 0 &&
			(c.spread != SPREAD_AGGRESSIVE || c.blocks > 1 || c.tf_change[i] < 0) {
			effectiveLowband = max(
				int(FREQ_BANDS[bandStart]),
				int(FREQ_BANDS[lowbandOffset]-FREQ_RANGE[i]),
			)
			fmt.Printf(
				"effectiveLowband %d off %d range %d\n",
				effectiveLowband, lowbandOffset, FREQ_RANGE[i],
			)

			foldstart := lowbandOffset
			for j := lowbandOffset - 1; j >= 0; j-- {
				if FREQ_BANDS[j] <= effectiveLowband {
					foldstart = j
					break
				}
			}

			foldend := lowbandOffset
			for j := lowbandOffset; j < MAX_BANDS; j++ {
				if FREQ_BANDS[j] >= effectiveLowband+FREQ_RANGE[i] {
					foldend = j
					break
				}
			}
			fmt.Printf("fold %d %d\n", foldstart, foldend)

			for j := foldstart; j < foldend; j++ {
				cm[0] |= int(c.frames[0].collapse_masks[j])
				cm[1] |= int(c.frames[c.stereo_pkt].collapse_masks[j])
			}
		} else {
			cm[0] = (1 << c.blocks) - 1
			cm[1] = cm[0]
		}

		fmt.Printf("cm %d %d\n", cm[0], cm[1])

		if c.dual_stereo && i == c.intensity_stereo {
			c.dual_stereo = false
			for j := FREQ_BANDS[bandStart] << lm; j < bandOffset; j++ {
				normMid[j] = (normMid[j] + normSide[j]) / 2.0
			}
		}

		var lowbandMid, lowbandMidOut []float32
		var lowbandSide, lowbandSideOut []float32

		if lowbandOffset != 0 &&
			(c.spread != SPREAD_AGGRESSIVE || c.blocks > 1 || c.tf_change[i] < 0) {
			lowbandMid = normMid[effectiveLowband<<lm : effectiveLowband<<lm+bandSize]
		}

		if i != bandEnd-1 {
			lowbandMidOut = normMid[bandOffset : bandOffset+bandSize]
		}

		if c.dual_stereo {
			if lowbandOffset != 0 &&
				(c.spread != SPREAD_AGGRESSIVE || c.blocks > 1 || c.tf_change[i] < 0) {
				lowbandSide = normSide[effectiveLowband<<lm : effectiveLowband<<lm+bandSize]
			}
			if i != bandEnd-1 {
				lowbandSideOut = normSide[bandOffset : bandOffset+bandSize]
			}

			cm[0] = c.DecodeBand(
				rd,
				i,
				x,
				nil,
				bandSize,
				b/2,
				c.blocks,
				lowbandMid,
				lowbandMidOut,
				lm,
				0,
				1.0,
				cm[0],
			)
			cm[1] = c.DecodeBand(
				rd,
				i,
				y,
				nil,
				bandSize,
				b/2,
				c.blocks,
				lowbandSide,
				lowbandSideOut,
				lm,
				0,
				1.0,
				cm[1],
			)
		} else {
			cm[0] = c.DecodeBand(
				rd,
				i,
				x,
				y,
				bandSize,
				b/2,
				c.blocks,
				lowbandMid,
				lowbandMidOut,
				lm,
				0,
				1.0,
				cm[0]|cm[1],
			)
			cm[1] = cm[0]
		}

		c.frames[0].collapse_masks[i] = uint8(cm[0])
		if c.stereo_pkt {
			c.frames[1].collapse_masks[i] = uint8(cm[1])
		}
		c.remaining += int32(c.pulses[i] + consumed)

		updateLowband = b > int(bandSize<<3)
	}
}

func (c *Celt) Decode(
	rd *comm.RangeDecoder,
	outBuf []float32,
	frameDuration comm.FrameDuration,
	bandStart, bandEnd int,
) {
	if bandEnd > MAX_BANDS {
		panic("bandEnd too large")
	}

	frameSize := int(frameDuration)

	// 计算 lm: log2(frameSize / SHORT_BLOCKSIZE) - 1
	ratio := frameSize / SHORT_BLOCKSIZE
	lm := 0
	for ratio > 1 {
		ratio >>= 1
		lm++
	}
	if lm > 0 {
		lm--
	}
	c.lm = lm

	fmt.Printf("framebits %d tell %d\n", rd.Len(), rd.Tell())

	silence := false
	if rd.Available() > 0 {
		silence = rd.DecodeLogP(15)
	} else {
		silence = true
	}

	fmt.Printf("silence %t\n", silence)

	if silence {
		rd.ToEnd() // 假设有 ToEnd 方法
	}

	c.ResetGains()

	if bandStart == 0 && rd.Available() >= 16 {
		c.ParsePostfilter(rd)
	}

	transient := false
	if c.lm != 0 && rd.Available() >= 3 {
		transient = rd.DecodeLogP(3)
	}

	fmt.Printf("duration %d, transient %t\n", c.lm, transient)

	if transient {
		c.blocks = 1 << c.lm
	} else {
		c.blocks = 1
	}
	c.blocksize = frameSize / c.blocks

	if !c.stereo_pkt {
		for i := range c.frames[0].energy {
			if c.frames[0].energy[i] < c.frames[1].energy[i] {
				c.frames[0].energy[i] = c.frames[1].energy[i]
			}
		}
	}

	for i := range c.frames {
		for j := range c.frames[i].collapse_masks {
			c.frames[i].collapse_masks[j] = 0
		}
	}

	c.DecodeCoarseEnergy(rd, bandStart, bandEnd)

	fmt.Printf("available %d tell %d frac %d\n", rd.Available(), rd.Tell(), rd.TellFrac())

	c.DecodeTfChanges(rd, bandStart, bandEnd, transient)

	fmt.Printf("available %d tell %d frac %d\n", rd.Available(), rd.Tell(), rd.TellFrac())

	c.DecodeAllocation(rd, bandStart, bandEnd)

	fmt.Printf("available %d tell %d frac %d\n", rd.Available(), rd.Tell(), rd.TellFrac())

	c.DecodeFineEnergy(rd, bandStart, bandEnd)

	coeff0 := make([]float32, MAX_FRAME_SIZE)
	coeff1 := make([]float32, MAX_FRAME_SIZE)

	fmt.Printf("available %d tell %d frac %d\n", rd.Available(), rd.Tell(), rd.TellFrac())

	c.DecodeBands(rd, bandStart, bandEnd, coeff0, coeff1)

	c.seed = uint32(rd.Range) // 假设 RangeDecoder 有 Range 字段
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
