package opsu

import (
	"github.com/project/silk/Inlines"
)

// BWExpander 带宽扩展器
type BWExpander struct{}

// Bwexpander32 对32位整数数组进行带宽扩展
func (b *BWExpander) Silk_bwexpander_32(
	ar []int32, // I/O AR滤波器系数 (无前导1)
	d int, // I 滤波器阶数
	chirp_Q16 int32, // I 带宽扩展因子 (Q16格式)
) {
	chirp_minus_one_Q16 := chirp_Q16 - 65536

	for i := 0; i < d-1; i++ {
		// 计算：ar[i] = (chirp_Q16 * ar[i]) >> 16
		ar[i] = Inlines.Silk_SMULWW(chirp_Q16, ar[i])

		// 更新chirp因子：chirp_Q16 += (chirp_Q16 * chirp_minus_one_Q16) >> 16
		chirp_Q16 += Inlines.Silk_RSHIFT_ROUND(
			Inlines.Silk_MUL(chirp_Q16, chirp_minus_one_Q16),
			16,
		)
	}
	ar[d-1] = Inlines.Silk_SMULWW(chirp_Q16, ar[d-1])
}

// Bwexpander 对16位整数数组进行带宽扩展
func (b *BWExpander) Silk_bwexpander(
	ar []int16, // I/O AR滤波器系数 (无前导1)
	d int, // I 滤波器阶数
	chirp_Q16 int32, // I 带宽扩展因子 (Q16格式)
) {
	chirp_minus_one_Q16 := chirp_Q16 - 65536

	for i := 0; i < d-1; i++ {
		// 计算：ar[i] = (chirp_Q16 * ar[i]) >> 16
		ar[i] = int16(Inlines.Silk_RSHIFT_ROUND(
			Inlines.Silk_MUL(chirp_Q16, int32(ar[i])),
			16,
		))

		// 更新chirp因子：chirp_Q16 += (chirp_Q16 * chirp_minus_one_Q16) >> 16
		chirp_Q16 += Inlines.Silk_RSHIFT_ROUND(
			Inlines.Silk_MUL(chirp_Q16, chirp_minus_one_Q16),
			16,
		)
	}
	ar[d-1] = int16(Inlines.Silk_RSHIFT_ROUND(
		Inlines.Silk_MUL(chirp_Q16, int32(ar[d-1])),
		16,
	))
}
