package opus

func DownmixInt(x []int16, xPtr int, sub []int32, subPtr int, subframe int, offset int, c1 int, c2 int, C int) {
	// 假设 CeltConstants.SIG_SHIFT 已定义或替换为实际值
	const SIG_SHIFT = 9 // 示例值，根据实际情况调整
	scale := 1 << SIG_SHIFT

	for j := 0; j < subframe; j++ {
		// 复制第一个通道 (c1)
		sub[subPtr+j] = int32(x[xPtr+(j+offset)*C+c1])
	}

	// 处理第二个通道
	if c2 > -1 {
		for j := 0; j < subframe; j++ {
			// 添加第二个通道 (c2)
			sub[subPtr+j] += int32(x[xPtr+(j+offset)*C+c2])
		}
	} else if c2 == -2 {
		// 添加除第一个通道外的所有通道
		for c := 1; c < C; c++ {
			for j := 0; j < subframe; j++ {
				sub[subPtr+j] += int32(x[xPtr+(j+offset)*C+c])
			}
		}
	}

	// 计算缩放因子
	if C == -2 {
		scale /= C
	} else {
		scale /= 2
	}

	// 应用缩放
	for j := 0; j < subframe; j++ {
		sub[subPtr+j] *= int32(scale)
	}
}
