package silk

import "github.com/dosgo/goOpus/comm"

// STAGE1 is the ICDF context for stage 1
var STAGE1 = &comm.ICDFContext{
	Total: 256,
	Dist: []int{
		7, 9, 10, 11, 12, 22, 46, 54, 55, 56, 59, 82, 174, 197, 200, 201, 202, 210, 234, 244, 245,
		246, 247, 249, 256,
	},
}

// STAGE2 is the ICDF context for stage 2
var STAGE2 = &comm.ICDFContext{
	Total: 256,
	Dist:  []int{85, 171, 256},
}

// STAGE3 is the ICDF context for stage 3
var STAGE3 = &comm.ICDFContext{
	Total: 256,
	Dist:  []int{51, 102, 154, 205, 256},
}

// MID_ONLY is the ICDF context for mid-only
var MID_ONLY = &comm.ICDFContext{
	Total: 256,
	Dist:  []int{192, 256},
}

// FRAME_TYPE_INACTIVE is the ICDF context for inactive frame type
var FRAME_TYPE_INACTIVE = &comm.ICDFContext{
	Total: 256,
	Dist:  []int{26, 256},
}

// FRAME_TYPE_ACTIVE is the ICDF context for active frame type
var FRAME_TYPE_ACTIVE = &comm.ICDFContext{
	Total: 256,
	Dist:  []int{24, 98, 246, 256},
}

// MSB_SUBFRAME_GAIN is the ICDF contexts for MSB subframe gain
var MSB_SUBFRAME_GAIN = []*comm.ICDFContext{
	{
		Total: 256,
		Dist:  []int{32, 144, 212, 241, 253, 254, 255, 256},
	},
	{
		Total: 256,
		Dist:  []int{2, 19, 64, 124, 186, 233, 252, 256},
	},
	{
		Total: 256,
		Dist:  []int{1, 4, 30, 101, 195, 245, 254, 256},
	},
}

// LSB_SUBFRAME_GAIN is the ICDF context for LSB subframe gain
var LSB_SUBFRAME_GAIN = &comm.ICDFContext{
	Total: 256,
	Dist:  []int{32, 64, 96, 128, 160, 192, 224, 256},
}

// DELTA_SUBFRAME_GAIN is the ICDF context for delta subframe gain
var DELTA_SUBFRAME_GAIN = &comm.ICDFContext{
	Total: 256,
	Dist: []int{
		6, 11, 22, 53, 185, 206, 214, 218, 221, 223, 225, 227, 228, 229, 230, 231, 232, 233, 234,
		235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252,
		253, 254, 255, 256,
	},
}

// LSF_STAGE1_NB_MB is the ICDF contexts for LSF stage 1 NB/MB
var LSF_STAGE1_NB_MB = []*comm.ICDFContext{
	{
		Total: 256,
		Dist: []int{
			44, 78, 108, 127, 148, 160, 171, 174, 177, 179, 195, 197, 199, 200, 205, 207, 208, 211,
			214, 215, 216, 218, 220, 222, 225, 226, 235, 244, 246, 253, 255, 256,
		},
	},
	{
		Total: 256,
		Dist: []int{
			1, 11, 12, 20, 23, 31, 39, 53, 66, 80, 81, 95, 107, 120, 131, 142, 154, 165, 175, 185,
			196, 204, 213, 221, 228, 236, 237, 238, 244, 245, 251, 256,
		},
	},
}

// LSF_STAGE1_WB is the ICDF contexts for LSF stage 1 WB
var LSF_STAGE1_WB = []*comm.ICDFContext{
	{
		Total: 256,
		Dist: []int{
			31, 52, 55, 72, 73, 81, 98, 102, 103, 121, 137, 141, 143, 146, 147, 157, 158, 161, 177,
			188, 204, 206, 208, 211, 213, 224, 225, 229, 238, 246, 253, 256,
		},
	},
	{
		Total: 256,
		Dist: []int{
			1, 5, 21, 26, 44, 55, 60, 74, 89, 90, 93, 105, 118, 132, 146, 152, 166, 178, 180, 186,
			187, 199, 211, 222, 232, 235, 245, 250, 251, 252, 253, 256,
		},
	},
}

type lsfStage2NbMbStruct struct {
	A   *comm.ICDFContext
	B   *comm.ICDFContext
	C   *comm.ICDFContext
	D   *comm.ICDFContext
	E   *comm.ICDFContext
	F   *comm.ICDFContext
	G   *comm.ICDFContext
	H   *comm.ICDFContext
	I   *comm.ICDFContext
	J   *comm.ICDFContext
	K   *comm.ICDFContext
	L   *comm.ICDFContext
	M   *comm.ICDFContext
	N   *comm.ICDFContext
	O   *comm.ICDFContext
	P   *comm.ICDFContext
	MAP [][]*comm.ICDFContext
}

var LSF_STAGE2_NB_MB = lsfStage2NbMbStruct{
	A: &comm.ICDFContext{
		Total: 256,
		Dist:  []int{1, 2, 3, 18, 242, 253, 254, 255, 256},
	},
	B: &comm.ICDFContext{
		Total: 256,
		Dist:  []int{1, 2, 4, 38, 221, 253, 254, 255, 256},
	},
	C: &comm.ICDFContext{
		Total: 256,
		Dist:  []int{1, 2, 6, 48, 197, 252, 254, 255, 256},
	},
	D: &comm.ICDFContext{
		Total: 256,
		Dist:  []int{1, 2, 10, 62, 185, 246, 254, 255, 256},
	},
	E: &comm.ICDFContext{
		Total: 256,
		Dist:  []int{1, 4, 20, 73, 174, 248, 254, 255, 256},
	},
	F: &comm.ICDFContext{
		Total: 256,
		Dist:  []int{1, 4, 21, 76, 166, 239, 254, 255, 256},
	},
	G: &comm.ICDFContext{
		Total: 256,
		Dist:  []int{1, 8, 32, 85, 159, 226, 252, 255, 256},
	},
	H: &comm.ICDFContext{
		Total: 256,
		Dist:  []int{1, 2, 20, 83, 161, 219, 249, 255, 256},
	},
}

// LSF_STAGE2_NB_MB contains the ICDF contexts for LSF stage 2 NB/MB
/*
var LSF_STAGE2_NB_MB = struct {
	A, B, C, D, E, F, G, H *comm.ICDFContext
	MAP                    [][]*comm.ICDFContext
}{
	A: &comm.ICDFContext{
		Total: 256,
		Dist:  []int{1, 2, 3, 18, 242, 253, 254, 255, 256},
	},
	B: &comm.ICDFContext{
		Total: 256,
		Dist:  []int{1, 2, 4, 38, 221, 253, 254, 255, 256},
	},
	C: &comm.ICDFContext{
		Total: 256,
		Dist:  []int{1, 2, 6, 48, 197, 252, 254, 255, 256},
	},
	D: &comm.ICDFContext{
		Total: 256,
		Dist:  []int{1, 2, 10, 62, 185, 246, 254, 255, 256},
	},
	E: &comm.ICDFContext{
		Total: 256,
		Dist:  []int{1, 4, 20, 73, 174, 248, 254, 255, 256},
	},
	F: &comm.ICDFContext{
		Total: 256,
		Dist:  []int{1, 4, 21, 76, 166, 239, 254, 255, 256},
	},
	G: &comm.ICDFContext{
		Total: 256,
		Dist:  []int{1, 8, 32, 85, 159, 226, 252, 255, 256},
	},
	H: &comm.ICDFContext{
		Total: 256,
		Dist:  []int{1, 2, 20, 83, 161, 219, 249, 255, 256},
	},
}
*/

func init() {
	LSF_STAGE2_NB_MB.MAP = [][]*comm.ICDFContext{
		{LSF_STAGE2_NB_MB.A, LSF_STAGE2_NB_MB.A, LSF_STAGE2_NB_MB.A, LSF_STAGE2_NB_MB.A, LSF_STAGE2_NB_MB.A, LSF_STAGE2_NB_MB.A, LSF_STAGE2_NB_MB.A, LSF_STAGE2_NB_MB.A, LSF_STAGE2_NB_MB.A, LSF_STAGE2_NB_MB.A}, // Placeholder for index 0
		{LSF_STAGE2_NB_MB.B, LSF_STAGE2_NB_MB.D, LSF_STAGE2_NB_MB.B, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.B, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.B, LSF_STAGE2_NB_MB.B, LSF_STAGE2_NB_MB.B},
		{LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.B, LSF_STAGE2_NB_MB.B, LSF_STAGE2_NB_MB.B, LSF_STAGE2_NB_MB.B, LSF_STAGE2_NB_MB.B, LSF_STAGE2_NB_MB.B, LSF_STAGE2_NB_MB.B, LSF_STAGE2_NB_MB.B, LSF_STAGE2_NB_MB.B},
		{LSF_STAGE2_NB_MB.B, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.B, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.B, LSF_STAGE2_NB_MB.B, LSF_STAGE2_NB_MB.B},
		{LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.D, LSF_STAGE2_NB_MB.D, LSF_STAGE2_NB_MB.D, LSF_STAGE2_NB_MB.D, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.C},
		{LSF_STAGE2_NB_MB.A, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.D, LSF_STAGE2_NB_MB.D, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.B, LSF_STAGE2_NB_MB.B},
		{LSF_STAGE2_NB_MB.A, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.B},
		{LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.D, LSF_STAGE2_NB_MB.G, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.F},
		{LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.G, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.E},
		{LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.H, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.E},
		{LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.D, LSF_STAGE2_NB_MB.D, LSF_STAGE2_NB_MB.D, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.D, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.C},
		{LSF_STAGE2_NB_MB.B, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.G, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.F},
		{LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.H, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.G, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.F},
		{LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.H, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.G, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.E},
		{LSF_STAGE2_NB_MB.D, LSF_STAGE2_NB_MB.D, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.E},
		{LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.D, LSF_STAGE2_NB_MB.D, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.E},
		{LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.G, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.F},
		{LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.G, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.E},
		{LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.H, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.F},
		{LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.G, LSF_STAGE2_NB_MB.H, LSF_STAGE2_NB_MB.G, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.G, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.E},
		{LSF_STAGE2_NB_MB.D, LSF_STAGE2_NB_MB.G, LSF_STAGE2_NB_MB.H, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.G, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.G, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.F},
		{LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.H, LSF_STAGE2_NB_MB.G, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.F},
		{LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.G, LSF_STAGE2_NB_MB.G, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.G, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.E},
		{LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.G, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.G, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.G, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.E},
		{LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.D, LSF_STAGE2_NB_MB.H, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.E},
		{LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.D, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.G, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.E},
		{LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.D, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.D, LSF_STAGE2_NB_MB.D, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.D, LSF_STAGE2_NB_MB.D, LSF_STAGE2_NB_MB.D},
		{LSF_STAGE2_NB_MB.B, LSF_STAGE2_NB_MB.B, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.D, LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.C},
		{LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.G, LSF_STAGE2_NB_MB.G, LSF_STAGE2_NB_MB.G, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.G, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.F},
		{LSF_STAGE2_NB_MB.D, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.D, LSF_STAGE2_NB_MB.D, LSF_STAGE2_NB_MB.C},
		{LSF_STAGE2_NB_MB.C, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.D, LSF_STAGE2_NB_MB.H, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.E},
		{LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.E, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.G, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.G, LSF_STAGE2_NB_MB.F, LSF_STAGE2_NB_MB.E},
	}
	lsf_stage2_wb.MAP = [][]*comm.ICDFContext{
		{lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I},
		{lsf_stage2_wb.K, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.K, lsf_stage2_wb.K, lsf_stage2_wb.K, lsf_stage2_wb.K, lsf_stage2_wb.K, lsf_stage2_wb.J, lsf_stage2_wb.J, lsf_stage2_wb.J, lsf_stage2_wb.I, lsf_stage2_wb.L},
		{lsf_stage2_wb.K, lsf_stage2_wb.N, lsf_stage2_wb.N, lsf_stage2_wb.L, lsf_stage2_wb.P, lsf_stage2_wb.M, lsf_stage2_wb.M, lsf_stage2_wb.N, lsf_stage2_wb.K, lsf_stage2_wb.N, lsf_stage2_wb.M, lsf_stage2_wb.N, lsf_stage2_wb.N, lsf_stage2_wb.M, lsf_stage2_wb.L, lsf_stage2_wb.L},
		{lsf_stage2_wb.I, lsf_stage2_wb.K, lsf_stage2_wb.J, lsf_stage2_wb.K, lsf_stage2_wb.K, lsf_stage2_wb.J, lsf_stage2_wb.J, lsf_stage2_wb.J, lsf_stage2_wb.J, lsf_stage2_wb.J, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.J},
		{lsf_stage2_wb.I, lsf_stage2_wb.O, lsf_stage2_wb.N, lsf_stage2_wb.M, lsf_stage2_wb.O, lsf_stage2_wb.M, lsf_stage2_wb.P, lsf_stage2_wb.N, lsf_stage2_wb.M, lsf_stage2_wb.M, lsf_stage2_wb.M, lsf_stage2_wb.N, lsf_stage2_wb.N, lsf_stage2_wb.M, lsf_stage2_wb.M, lsf_stage2_wb.L},
		{lsf_stage2_wb.I, lsf_stage2_wb.L, lsf_stage2_wb.N, lsf_stage2_wb.N, lsf_stage2_wb.M, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.N, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.K, lsf_stage2_wb.M},
		{lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I},
		{lsf_stage2_wb.I, lsf_stage2_wb.K, lsf_stage2_wb.O, lsf_stage2_wb.L, lsf_stage2_wb.P, lsf_stage2_wb.K, lsf_stage2_wb.N, lsf_stage2_wb.L, lsf_stage2_wb.M, lsf_stage2_wb.N, lsf_stage2_wb.N, lsf_stage2_wb.M, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.K, lsf_stage2_wb.L},
		{lsf_stage2_wb.I, lsf_stage2_wb.O, lsf_stage2_wb.K, lsf_stage2_wb.O, lsf_stage2_wb.O, lsf_stage2_wb.M, lsf_stage2_wb.N, lsf_stage2_wb.M, lsf_stage2_wb.O, lsf_stage2_wb.N, lsf_stage2_wb.M, lsf_stage2_wb.M, lsf_stage2_wb.N, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.L},
		{lsf_stage2_wb.K, lsf_stage2_wb.J, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I},
		{lsf_stage2_wb.I, lsf_stage2_wb.J, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.J},
		{lsf_stage2_wb.K, lsf_stage2_wb.K, lsf_stage2_wb.L, lsf_stage2_wb.M, lsf_stage2_wb.N, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.K, lsf_stage2_wb.K, lsf_stage2_wb.J, lsf_stage2_wb.L},
		{lsf_stage2_wb.K, lsf_stage2_wb.K, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.M, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.K, lsf_stage2_wb.J, lsf_stage2_wb.L},
		{lsf_stage2_wb.L, lsf_stage2_wb.M, lsf_stage2_wb.M, lsf_stage2_wb.M, lsf_stage2_wb.O, lsf_stage2_wb.M, lsf_stage2_wb.M, lsf_stage2_wb.N, lsf_stage2_wb.L, lsf_stage2_wb.N, lsf_stage2_wb.M, lsf_stage2_wb.M, lsf_stage2_wb.N, lsf_stage2_wb.M, lsf_stage2_wb.L, lsf_stage2_wb.M},
		{lsf_stage2_wb.I, lsf_stage2_wb.O, lsf_stage2_wb.M, lsf_stage2_wb.N, lsf_stage2_wb.M, lsf_stage2_wb.P, lsf_stage2_wb.N, lsf_stage2_wb.K, lsf_stage2_wb.O, lsf_stage2_wb.N, lsf_stage2_wb.P, lsf_stage2_wb.M, lsf_stage2_wb.M, lsf_stage2_wb.L, lsf_stage2_wb.N, lsf_stage2_wb.L},
		{lsf_stage2_wb.I, lsf_stage2_wb.J, lsf_stage2_wb.I, lsf_stage2_wb.J, lsf_stage2_wb.J, lsf_stage2_wb.J, lsf_stage2_wb.J, lsf_stage2_wb.J, lsf_stage2_wb.J, lsf_stage2_wb.J, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.J, lsf_stage2_wb.I},
		{lsf_stage2_wb.J, lsf_stage2_wb.O, lsf_stage2_wb.N, lsf_stage2_wb.P, lsf_stage2_wb.N, lsf_stage2_wb.M, lsf_stage2_wb.N, lsf_stage2_wb.L, lsf_stage2_wb.M, lsf_stage2_wb.N, lsf_stage2_wb.M, lsf_stage2_wb.M, lsf_stage2_wb.M, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.M},
		{lsf_stage2_wb.J, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.M, lsf_stage2_wb.M, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.N, lsf_stage2_wb.K, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.N, lsf_stage2_wb.N, lsf_stage2_wb.N, lsf_stage2_wb.L, lsf_stage2_wb.M},
		{lsf_stage2_wb.K, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.K, lsf_stage2_wb.K, lsf_stage2_wb.K, lsf_stage2_wb.L, lsf_stage2_wb.K, lsf_stage2_wb.J, lsf_stage2_wb.K, lsf_stage2_wb.J, lsf_stage2_wb.K, lsf_stage2_wb.J, lsf_stage2_wb.J, lsf_stage2_wb.J, lsf_stage2_wb.M},
		{lsf_stage2_wb.I, lsf_stage2_wb.K, lsf_stage2_wb.L, lsf_stage2_wb.N, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.K, lsf_stage2_wb.K, lsf_stage2_wb.K, lsf_stage2_wb.J, lsf_stage2_wb.J, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I},
		{lsf_stage2_wb.L, lsf_stage2_wb.M, lsf_stage2_wb.L, lsf_stage2_wb.N, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.K, lsf_stage2_wb.K, lsf_stage2_wb.J, lsf_stage2_wb.J, lsf_stage2_wb.J, lsf_stage2_wb.J, lsf_stage2_wb.J, lsf_stage2_wb.K, lsf_stage2_wb.K, lsf_stage2_wb.M},
		{lsf_stage2_wb.K, lsf_stage2_wb.O, lsf_stage2_wb.L, lsf_stage2_wb.P, lsf_stage2_wb.P, lsf_stage2_wb.M, lsf_stage2_wb.N, lsf_stage2_wb.M, lsf_stage2_wb.N, lsf_stage2_wb.L, lsf_stage2_wb.N, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.K, lsf_stage2_wb.L, lsf_stage2_wb.L},
		{lsf_stage2_wb.K, lsf_stage2_wb.L, lsf_stage2_wb.N, lsf_stage2_wb.O, lsf_stage2_wb.O, lsf_stage2_wb.L, lsf_stage2_wb.N, lsf_stage2_wb.L, lsf_stage2_wb.M, lsf_stage2_wb.M, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.K, lsf_stage2_wb.M},
		{lsf_stage2_wb.J, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.M, lsf_stage2_wb.M, lsf_stage2_wb.M, lsf_stage2_wb.M, lsf_stage2_wb.L, lsf_stage2_wb.N, lsf_stage2_wb.N, lsf_stage2_wb.N, lsf_stage2_wb.L, lsf_stage2_wb.J, lsf_stage2_wb.J, lsf_stage2_wb.J, lsf_stage2_wb.J},
		{lsf_stage2_wb.K, lsf_stage2_wb.N, lsf_stage2_wb.L, lsf_stage2_wb.O, lsf_stage2_wb.O, lsf_stage2_wb.M, lsf_stage2_wb.P, lsf_stage2_wb.M, lsf_stage2_wb.M, lsf_stage2_wb.N, lsf_stage2_wb.L, lsf_stage2_wb.M, lsf_stage2_wb.M, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.L},
		{lsf_stage2_wb.I, lsf_stage2_wb.O, lsf_stage2_wb.J, lsf_stage2_wb.J, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I},
		{lsf_stage2_wb.I, lsf_stage2_wb.O, lsf_stage2_wb.O, lsf_stage2_wb.L, lsf_stage2_wb.N, lsf_stage2_wb.K, lsf_stage2_wb.N, lsf_stage2_wb.N, lsf_stage2_wb.L, lsf_stage2_wb.M, lsf_stage2_wb.M, lsf_stage2_wb.P, lsf_stage2_wb.P, lsf_stage2_wb.M, lsf_stage2_wb.M, lsf_stage2_wb.M},
		{lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.P, lsf_stage2_wb.L, lsf_stage2_wb.N, lsf_stage2_wb.M, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.K, lsf_stage2_wb.K, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.K, lsf_stage2_wb.L},
		{lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.J, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.K, lsf_stage2_wb.J, lsf_stage2_wb.K, lsf_stage2_wb.J, lsf_stage2_wb.J, lsf_stage2_wb.K, lsf_stage2_wb.K, lsf_stage2_wb.K, lsf_stage2_wb.J, lsf_stage2_wb.J},
		{lsf_stage2_wb.I, lsf_stage2_wb.L, lsf_stage2_wb.K, lsf_stage2_wb.N, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.K, lsf_stage2_wb.L, lsf_stage2_wb.K, lsf_stage2_wb.J, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.J, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.J},
		{lsf_stage2_wb.L, lsf_stage2_wb.N, lsf_stage2_wb.N, lsf_stage2_wb.M, lsf_stage2_wb.P, lsf_stage2_wb.N, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.K, lsf_stage2_wb.L, lsf_stage2_wb.K, lsf_stage2_wb.K, lsf_stage2_wb.J, lsf_stage2_wb.I, lsf_stage2_wb.J, lsf_stage2_wb.I},
		{lsf_stage2_wb.K, lsf_stage2_wb.L, lsf_stage2_wb.N, lsf_stage2_wb.L, lsf_stage2_wb.M, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.L, lsf_stage2_wb.K, lsf_stage2_wb.J, lsf_stage2_wb.K, lsf_stage2_wb.O, lsf_stage2_wb.M, lsf_stage2_wb.I, lsf_stage2_wb.I, lsf_stage2_wb.I},
	}
}

var lsf_stage2_wb = lsfStage2NbMbStruct{
	// 创建 ICDFContext 实例
	I: &comm.ICDFContext{
		Total: 256,
		Dist:  []int{143, 193, 256},
	},
	J: &comm.ICDFContext{
		Total: 256,
		Dist:  []int{68, 80, 101, 118, 137, 159, 189, 213, 230, 246, 256},
	},
	K: &comm.ICDFContext{
		Total: 256,
		Dist:  []int{91, 137, 176, 195, 209, 221, 229, 236, 242, 247, 252, 256},
	},
	L: &comm.ICDFContext{
		Total: 256,
		Dist:  []int{33, 55, 73, 89, 104, 118, 132, 145, 158, 168, 177, 186, 194, 200, 206, 212, 217, 221, 225, 229, 232, 235, 238, 240, 242, 244, 246, 248, 250, 252, 253, 254, 255, 256},
	},
	M: &comm.ICDFContext{
		Total: 256,
		Dist:  []int{77, 157, 256},
	},
	N: &comm.ICDFContext{
		Total: 256,
		Dist:  []int{185, 200, 213, 226, 235, 244, 250, 256},
	},
	O: &comm.ICDFContext{
		Total: 256,
		Dist:  []int{57, 91, 112, 132, 147, 160, 172, 185, 195, 205, 214, 224, 233, 241, 248, 256},
	},
	P: &comm.ICDFContext{
		Total: 256,
		Dist:  []int{15, 31, 45, 57, 69, 81, 92, 103, 114, 124, 133, 142, 151, 160, 168, 176, 184, 192, 199, 206, 212, 218, 223, 227, 232, 236, 240, 244, 247, 251, 254, 256},
	},
}

// LSF_STAGE2_EXTENSION is the ICDF context for LSF stage 2 extension
var LSF_STAGE2_EXTENSION = &comm.ICDFContext{
	Total: 256,
	Dist:  []int{156, 216, 240, 249, 253, 255, 256},
}

// LSF_PRED_WEIGHT_NB_MB contains the prediction weights for NB/MB
var LSF_PRED_WEIGHT_NB_MB = [][]uint8{
	{179, 138, 140, 148, 151, 149, 153, 151, 163},
	{116, 67, 82, 59, 92, 72, 100, 89, 92},
}

// LSF_PRED_WEIGHT_WB contains the prediction weights for WB
var LSF_PRED_WEIGHT_WB = [][]uint8{
	{175, 148, 160, 176, 178, 173, 174, 164, 177, 174, 196, 182, 198, 192, 182},
	{68, 62, 66, 60, 72, 117, 85, 90, 118, 136, 151, 142, 160, 142, 155},
}

// LSF_PRED_WEIGHT_INDEX_NB_MB contains the prediction weight indices for NB/MB
var LSF_PRED_WEIGHT_INDEX_NB_MB = [][]int{
	{0, 1, 0, 0, 0, 0, 0, 0, 0},
	{1, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0},
	{1, 1, 1, 0, 0, 0, 0, 1, 0},
	{0, 1, 0, 0, 0, 0, 0, 0, 0},
	{0, 1, 0, 0, 0, 0, 0, 0, 0},
	{1, 0, 1, 1, 0, 0, 0, 1, 0},
	{0, 1, 1, 0, 0, 1, 1, 0, 0},
	{0, 0, 1, 1, 0, 1, 0, 1, 1},
	{0, 0, 1, 1, 0, 0, 1, 1, 1},
	{0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 1, 0, 1, 1, 1, 1, 1, 0},
	{0, 1, 0, 1, 1, 1, 1, 1, 0},
	{0, 1, 1, 1, 1, 1, 1, 1, 0},
	{1, 0, 1, 1, 0, 1, 1, 1, 1},
	{0, 1, 1, 1, 1, 1, 0, 1, 0},
	{0, 0, 1, 1, 0, 1, 0, 1, 0},
	{0, 0, 1, 1, 1, 0, 1, 1, 1},
	{0, 1, 1, 0, 0, 1, 1, 1, 0},
	{0, 0, 0, 1, 1, 1, 0, 1, 0},
	{0, 1, 1, 0, 0, 1, 0, 1, 0},
	{0, 1, 1, 0, 0, 0, 1, 1, 0},
	{0, 0, 0, 0, 0, 1, 1, 1, 1},
	{0, 0, 1, 1, 0, 0, 0, 1, 1},
	{0, 0, 0, 1, 0, 1, 1, 1, 1},
	{0, 1, 1, 1, 1, 1, 1, 1, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 1, 0, 1, 1, 0, 1, 0},
	{1, 0, 0, 1, 0, 0, 0, 0, 0},
	{0, 0, 0, 1, 1, 0, 1, 0, 1},
	{1, 0, 1, 1, 0, 1, 1, 1, 1},
}

// LSF_PRED_WEIGHT_INDEX_WB contains the prediction weight indices for WB
var LSF_PRED_WEIGHT_INDEX_WB = [][]int{
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
	{0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0},
	{0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0},
	{0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1},
	{0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
	{0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0},
	{0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0},
	{0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0},
	{0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1},
	{0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0},
	{0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0},
	{0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0},
	{0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0},
	{0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0},
	{0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
	{0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1},
	{0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
	{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
	{0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0},
}

// LSF_CODEBOOK_NB_MB contains the LSF codebook for NB/MB
var LSF_CODEBOOK_NB_MB = [][]uint8{
	{12, 35, 60, 83, 108, 132, 157, 180, 206, 228},
	{15, 32, 55, 77, 101, 125, 151, 175, 201, 225},
	{19, 42, 66, 89, 114, 137, 162, 184, 209, 230},
	{12, 25, 50, 72, 97, 120, 147, 172, 200, 223},
	{26, 44, 69, 90, 114, 135, 159, 180, 205, 225},
	{13, 22, 53, 80, 106, 130, 156, 180, 205, 228},
	{15, 25, 44, 64, 90, 115, 142, 168, 196, 222},
	{19, 24, 62, 82, 100, 120, 145, 168, 190, 214},
	{22, 31, 50, 79, 103, 120, 151, 170, 203, 227},
	{21, 29, 45, 65, 106, 124, 150, 171, 196, 224},
	{30, 49, 75, 97, 121, 142, 165, 186, 209, 229},
	{19, 25, 52, 70, 93, 116, 143, 166, 192, 219},
	{26, 34, 62, 75, 97, 118, 145, 167, 194, 217},
	{25, 33, 56, 70, 91, 113, 143, 165, 196, 223},
	{21, 34, 51, 72, 97, 117, 145, 171, 196, 222},
	{20, 29, 50, 67, 90, 117, 144, 168, 197, 221},
	{22, 31, 48, 66, 95, 117, 146, 168, 196, 222},
	{24, 33, 51, 77, 116, 134, 158, 180, 200, 224},
	{21, 28, 70, 87, 106, 124, 149, 170, 194, 217},
	{26, 33, 53, 64, 83, 117, 152, 173, 204, 225},
	{27, 34, 65, 95, 108, 129, 155, 174, 210, 225},
	{20, 26, 72, 99, 113, 131, 154, 176, 200, 219},
	{34, 43, 61, 78, 93, 114, 155, 177, 205, 229},
	{23, 29, 54, 97, 124, 138, 163, 179, 209, 229},
	{30, 38, 56, 89, 118, 129, 158, 178, 200, 231},
	{21, 29, 49, 63, 85, 111, 142, 163, 193, 222},
	{27, 48, 77, 103, 133, 158, 179, 196, 215, 232},
	{29, 47, 74, 99, 124, 151, 176, 198, 220, 237},
	{33, 42, 61, 76, 93, 121, 155, 174, 207, 225},
	{29, 53, 87, 112, 136, 154, 170, 188, 208, 227},
	{24, 30, 52, 84, 131, 150, 166, 186, 203, 229},
	{37, 48, 64, 84, 104, 118, 156, 177, 201, 230},
}

// LSF_CODEBOOK_WB contains the LSF codebook for WB
var LSF_CODEBOOK_WB = [][]uint8{
	{7, 23, 38, 54, 69, 85, 100, 116, 131, 147, 162, 178, 193, 208, 223, 239},
	{13, 25, 41, 55, 69, 83, 98, 112, 127, 142, 157, 171, 187, 203, 220, 236},
	{15, 21, 34, 51, 61, 78, 92, 106, 126, 136, 152, 167, 185, 205, 225, 240},
	{10, 21, 36, 50, 63, 79, 95, 110, 126, 141, 157, 173, 189, 205, 221, 237},
	{17, 20, 37, 51, 59, 78, 89, 107, 123, 134, 150, 164, 184, 205, 224, 240},
	{10, 15, 32, 51, 67, 81, 96, 112, 129, 142, 158, 173, 189, 204, 220, 236},
	{8, 21, 37, 51, 65, 79, 98, 113, 126, 138, 155, 168, 179, 192, 209, 218},
	{12, 15, 34, 55, 63, 78, 87, 108, 118, 131, 148, 167, 185, 203, 219, 236},
	{16, 19, 32, 36, 56, 79, 91, 108, 118, 136, 154, 171, 186, 204, 220, 237},
	{11, 28, 43, 58, 74, 89, 105, 120, 135, 150, 165, 180, 196, 211, 226, 241},
	{6, 16, 33, 46, 60, 75, 92, 107, 123, 137, 156, 169, 185, 199, 214, 225},
	{11, 19, 30, 44, 57, 74, 89, 105, 121, 135, 152, 169, 186, 202, 218, 234},
	{12, 19, 29, 46, 57, 71, 88, 100, 120, 132, 148, 165, 182, 199, 216, 233},
	{17, 23, 35, 46, 56, 77, 92, 106, 123, 134, 152, 167, 185, 204, 222, 237},
	{14, 17, 45, 53, 63, 75, 89, 107, 115, 132, 151, 171, 188, 206, 221, 240},
	{9, 16, 29, 40, 56, 71, 88, 103, 119, 137, 154, 171, 189, 205, 222, 237},
	{16, 19, 36, 48, 57, 76, 87, 105, 118, 132, 150, 167, 185, 202, 218, 236},
	{12, 17, 29, 54, 71, 81, 94, 104, 126, 136, 149, 164, 182, 201, 221, 237},
	{15, 28, 47, 62, 79, 97, 115, 129, 142, 155, 168, 180, 194, 208, 223, 238},
	{8, 14, 30, 45, 62, 78, 94, 111, 127, 143, 159, 175, 192, 207, 223, 239},
	{17, 30, 49, 62, 79, 92, 107, 119, 132, 145, 160, 174, 190, 204, 220, 235},
	{14, 19, 36, 45, 61, 76, 91, 108, 121, 138, 154, 172, 189, 205, 222, 238},
	{12, 18, 31, 45, 60, 76, 91, 107, 123, 138, 154, 171, 187, 204, 221, 236},
	{13, 17, 31, 43, 53, 70, 83, 103, 114, 131, 149, 167, 185, 203, 220, 237},
	{17, 22, 35, 42, 58, 78, 93, 110, 125, 139, 155, 170, 188, 206, 224, 240},
	{8, 15, 34, 50, 67, 83, 99, 115, 131, 146, 162, 178, 193, 209, 224, 239},
	{13, 16, 41, 66, 73, 86, 95, 111, 128, 137, 150, 163, 183, 206, 225, 241},
	{17, 25, 37, 52, 63, 75, 92, 102, 119, 132, 144, 160, 175, 191, 212, 231},
	{19, 31, 49, 65, 83, 100, 117, 133, 147, 161, 174, 187, 200, 213, 227, 242},
	{18, 31, 52, 68, 88, 103, 117, 126, 138, 149, 163, 177, 192, 207, 223, 239},
	{16, 29, 47, 61, 76, 90, 106, 119, 133, 147, 161, 176, 193, 209, 224, 240},
	{15, 21, 35, 50, 61, 73, 86, 97, 110, 119, 129, 141, 175, 198, 218, 237},
}

// LSF_WEIGHT_NB_MB contains the LSF weights for NB/MB
var LSF_WEIGHT_NB_MB = [][]uint16{
	{2897, 2314, 2314, 2314, 2287, 2287, 2314, 2300, 2327, 2287},
	{2888, 2580, 2394, 2367, 2314, 2274, 2274, 2274, 2274, 2194},
	{2487, 2340, 2340, 2314, 2314, 2314, 2340, 2340, 2367, 2354},
	{3216, 2766, 2340, 2340, 2314, 2274, 2221, 2207, 2261, 2194},
	{2460, 2474, 2367, 2394, 2394, 2394, 2394, 2367, 2407, 2314},
	{3479, 3056, 2127, 2207, 2274, 2274, 2274, 2287, 2314, 2261},
	{3282, 3141, 2580, 2394, 2247, 2221, 2207, 2194, 2194, 2114},
	{4096, 3845, 2221, 2620, 2620, 2407, 2314, 2394, 2367, 2074},
	{3178, 3244, 2367, 2221, 2553, 2434, 2340, 2314, 2167, 2221},
	{3338, 3488, 2726, 2194, 2261, 2460, 2354, 2367, 2207, 2101},
	{2354, 2420, 2327, 2367, 2394, 2420, 2420, 2420, 2460, 2367},
	{3779, 3629, 2434, 2527, 2367, 2274, 2274, 2300, 2207, 2048},
	{3254, 3225, 2713, 2846, 2447, 2327, 2300, 2300, 2274, 2127},
	{3263, 3300, 2753, 2806, 2447, 2261, 2261, 2247, 2127, 2101},
	{2873, 2981, 2633, 2367, 2407, 2354, 2194, 2247, 2247, 2114},
	{3225, 3197, 2633, 2580, 2274, 2181, 2247, 2221, 2221, 2141},
	{3178, 3310, 2740, 2407, 2274, 2274, 2274, 2287, 2194, 2114},
	{3141, 3272, 2460, 2061, 2287, 2500, 2367, 2487, 2434, 2181},
	{3507, 3282, 2314, 2700, 2647, 2474, 2367, 2394, 2340, 2127},
	{3423, 3535, 3038, 3056, 2300, 1950, 2221, 2274, 2274, 2274},
	{3404, 3366, 2087, 2687, 2873, 2354, 2420, 2274, 2474, 2540},
	{3760, 3488, 1950, 2660, 2897, 2527, 2394, 2367, 2460, 2261},
	{3028, 3272, 2740, 2888, 2740, 2154, 2127, 2287, 2234, 2247},
	{3695, 3657, 2025, 1969, 2660, 2700, 2580, 2500, 2327, 2367},
	{3207, 3413, 2354, 2074, 2888, 2888, 2340, 2487, 2247, 2167},
	{3338, 3366, 2846, 2780, 2327, 2154, 2274, 2287, 2114, 2061},
	{2327, 2300, 2181, 2167, 2181, 2367, 2633, 2700, 2700, 2553},
	{2407, 2434, 2221, 2261, 2221, 2221, 2340, 2420, 2607, 2700},
	{3038, 3244, 2806, 2888, 2474, 2074, 2300, 2314, 2354, 2380},
	{2221, 2154, 2127, 2287, 2500, 2793, 2793, 2620, 2580, 2367},
	{3676, 3713, 2234, 1838, 2181, 2753, 2726, 2673, 2513, 2207},
	{2793, 3160, 2726, 2553, 2846, 2513, 2181, 2394, 2221, 2181},
}

// LSF_WEIGHT_WB contains the LSF weights for WB
var LSF_WEIGHT_WB = [][]uint16{
	{3657, 2925, 2925, 2925, 2925, 2925, 2925, 2925, 2925, 2925, 2925, 2925, 2963, 2963, 2925, 2846},
	{3216, 3085, 2972, 3056, 3056, 3010, 3010, 3010, 2963, 2963, 3010, 2972, 2888, 2846, 2846, 2726},
	{3920, 4014, 2981, 3207, 3207, 2934, 3056, 2846, 3122, 3244, 2925, 2846, 2620, 2553, 2780, 2925},
	{3516, 3197, 3010, 3103, 3019, 2888, 2925, 2925, 2925, 2925, 2888, 2888, 2888, 2888, 2888, 2753},
	{5054, 5054, 2934, 3573, 3385, 3056, 3085, 2793, 3160, 3160, 2972, 2846, 2513, 2540, 2753, 2888},
	{4428, 4149, 2700, 2753, 2972, 3010, 2925, 2846, 2981, 3019, 2925, 2925, 2925, 2925, 2888, 2726},
	{3620, 3019, 2972, 3056, 3056, 2873, 2806, 3056, 3216, 3047, 2981, 3291, 3291, 2981, 3310, 2991},
	{5227, 5014, 2540, 3338, 3526, 3385, 3197, 3094, 3376, 2981, 2700, 2647, 2687, 2793, 2846, 2673},
	{5081, 5174, 4615, 4428, 2460, 2897, 3047, 3207, 3169, 2687, 2740, 2888, 2846, 2793, 2846, 2700},
	{3122, 2888, 2963, 2925, 2925, 2925, 2925, 2963, 2963, 2963, 2963, 2925, 2925, 2963, 2963, 2963},
	{4202, 3207, 2981, 3103, 3010, 2888, 2888, 2925, 2972, 2873, 2916, 3019, 2972, 3010, 3197, 2873},
	{3760, 3760, 3244, 3103, 2981, 2888, 2925, 2888, 2972, 2934, 2793, 2793, 2846, 2888, 2888, 2660},
	{3854, 4014, 3207, 3122, 3244, 2934, 3047, 2963, 2963, 3085, 2846, 2793, 2793, 2793, 2793, 2580},
	{3845, 4080, 3357, 3516, 3094, 2740, 3010, 2934, 3122, 3085, 2846, 2846, 2647, 2647, 2846, 2806},
	{5147, 4894, 3225, 3845, 3441, 3169, 2897, 3413, 3451, 2700, 2580, 2673, 2740, 2846, 2806, 2753},
	{4109, 3789, 3291, 3160, 2925, 2888, 2888, 2925, 2793, 2740, 2793, 2740, 2793, 2846, 2888, 2806},
	{5081, 5054, 3047, 3545, 3244, 3056, 3085, 2944, 3103, 2897, 2740, 2740, 2740, 2846, 2793, 2620},
	{4309, 4309, 2860, 2527, 3207, 3376, 3376, 3075, 3075, 3376, 3056, 2846, 2647, 2580, 2726, 2753},
	{3056, 2916, 2806, 2888, 2740, 2687, 2897, 3103, 3150, 3150, 3216, 3169, 3056, 3010, 2963, 2846},
	{4375, 3882, 2925, 2888, 2846, 2888, 2846, 2846, 2888, 2888, 2888, 2846, 2888, 2925, 2888, 2846},
	{2981, 2916, 2916, 2981, 2981, 3056, 3122, 3216, 3150, 3056, 3010, 2972, 2972, 2972, 2925, 2740},
	{4229, 4149, 3310, 3347, 2925, 2963, 2888, 2981, 2981, 2846, 2793, 2740, 2846, 2846, 2846, 2793},
	{4080, 4014, 3103, 3010, 2925, 2925, 2925, 2888, 2925, 2925, 2846, 2846, 2846, 2793, 2888, 2780},
	{4615, 4575, 3169, 3441, 3207, 2981, 2897, 3038, 3122, 2740, 2687, 2687, 2687, 2740, 2793, 2700},
	{4149, 4269, 3789, 3657, 2726, 2780, 2888, 2888, 3010, 2972, 2925, 2846, 2687, 2687, 2793, 2888},
	{4215, 3554, 2753, 2846, 2846, 2888, 2888, 2888, 2925, 2925, 2888, 2925, 2925, 2925, 2963, 2888},
	{5174, 4921, 2261, 3432, 3789, 3479, 3347, 2846, 3310, 3479, 3150, 2897, 2460, 2487, 2753, 2925},
	{3451, 3685, 3122, 3197, 3357, 3047, 3207, 3207, 2981, 3216, 3085, 2925, 2925, 2687, 2540, 2434},
	{2981, 3010, 2793, 2793, 2740, 2793, 2846, 2972, 3056, 3103, 3150, 3150, 3150, 3103, 3010, 3010},
	{2944, 2873, 2687, 2726, 2780, 3010, 3432, 3545, 3357, 3244, 3056, 3010, 2963, 2925, 2888, 2846},
	{3019, 2944, 2897, 3010, 3010, 2972, 3019, 3103, 3056, 3056, 3010, 2888, 2846, 2925, 2925, 2888},
	{3920, 3967, 3010, 3197, 3357, 3216, 3291, 3291, 3479, 3704, 3441, 2726, 2181, 2460, 2580, 2607},
}

// LSF_MIN_SPACING_NB_MB contains the minimum spacing for NB/MB
var LSF_MIN_SPACING_NB_MB = []int16{250, 3, 6, 3, 3, 3, 4, 3, 3, 3, 461}

// LSF_MIN_SPACING_WB contains the minimum spacing for WB
var LSF_MIN_SPACING_WB = []int16{100, 3, 40, 3, 3, 3, 5, 14, 14, 10, 11, 3, 8, 9, 7, 3, 347}

// LSF_INTERPOLATION_INDEX is the ICDF context for LSF interpolation index
var LSF_INTERPOLATION_INDEX = &comm.ICDFContext{
	Total: 256,
	Dist:  []int{13, 35, 64, 75, 256},
}

// LSF_ORDERING_NB_MB contains the ordering for NB/MB
var LSF_ORDERING_NB_MB = []uint8{0, 9, 6, 3, 4, 5, 8, 1, 2, 7}

// LSF_ORDERING_WB contains the ordering for WB
var LSF_ORDERING_WB = []uint8{0, 15, 8, 7, 4, 11, 12, 3, 2, 13, 10, 5, 6, 9, 14, 1}

// COSINE contains cosine values
var COSINE = []int16{
	4096, 4095, 4091, 4085, 4076, 4065, 4052, 4036, 4017, 3997, 3973, 3948, 3920, 3889, 3857, 3822,
	3784, 3745, 3703, 3659, 3613, 3564, 3513, 3461, 3406, 3349, 3290, 3229, 3166, 3102, 3035, 2967,
	2896, 2824, 2751, 2676, 2599, 2520, 2440, 2359, 2276, 2191, 2106, 2019, 1931, 1842, 1751, 1660,
	1568, 1474, 1380, 1285, 1189, 1093, 995, 897, 799, 700, 601, 501, 401, 301, 201, 101, 0, -101,
	-201, -301, -401, -501, -601, -700, -799, -897, -995, -1093, -1189, -1285, -1380, -1474, -1568,
	-1660, -1751, -1842, -1931, -2019, -2106, -2191, -2276, -2359, -2440, -2520, -2599, -2676,
	-2751, -2824, -2896, -2967, -3035, -3102, -3166, -3229, -3290, -3349, -3406, -3461, -3513,
	-3564, -3613, -3659, -3703, -3745, -3784, -3822, -3857, -3889, -3920, -3948, -3973, -3997,
	-4017, -4036, -4052, -4065, -4076, -4085, -4091, -4095, -4096,
}

// PITCH_DELTA is the ICDF context for pitch delta
var PITCH_DELTA = &comm.ICDFContext{
	Total: 256,
	Dist: []int{
		46, 48, 50, 53, 57, 63, 73, 88, 114, 152, 182, 204, 219, 229, 236, 242, 246, 250, 252, 254,
		256,
	},
}
