package kissfft

const (
	MAXFACTORS    = 8
	SIG_SHIFT     = 15
	Q15_SCALE     = 1 << 15
	SQRT_HALF_Q15 = 23170 // Q15(0.7071067812)
	Q15_28378     = -28378
)

// FFTState 存储FFT配置状态
type FFTState struct {
	Factors    [2 * MAXFACTORS]int // Radix factors (p and m)
	Bitrev     []int               // Bit-reversal permutation
	Twiddles   []int               // Complex twiddle factors
	Nfft       int                 // FFT size
	Scale      int                 // Scaling value
	ScaleShift int                 // Scaling shift count
	Shift      int                 // Overall shift
}

// NewFFTState 创建指定大小的FFT状态
func NewFFTState(nfft int, inverse bool) *FFTState {
	// 简化实现: 实际应用中需填充factors和bitrev
	return &FFTState{
		Nfft:       nfft,
		Scale:      Q15_SCALE,
		ScaleShift: 1,
		Twiddles:   make([]int, 2*nfft),
		Bitrev:     make([]int, nfft),
	}
}

// MULT16_32_Q15 Q15定点数乘法
func MULT16_32_Q15(a int, b int) int {
	return (a * b) >> SIG_SHIFT
}

// SHR32 带符号右移
func SHR32(a, shift int) int {
	if shift > 0 {
		return a >> shift
	}
	return a << (-shift)
}

// sMul Q15乘法包装函数
func sMul(a, b int) int {
	return MULT16_32_Q15(b, a)
}

// halfOf 折半计算
func halfOf(x int) int {
	return x >> 1
}

// kf_bfly2 基2蝴蝶运算
func kf_bfly2(fout []int, fout_ptr, m, N int) {
	fout_ptr_base := fout_ptr
	tw := SQRT_HALF_Q15 // Q15(0.7071067812)

	for i := 0; i < N; i++ {
		fout_ptr = fout_ptr_base + i*8
		Fout2 := fout_ptr + 8

		// 第一组蝶形
		t_r := fout[Fout2+0]
		t_i := fout[Fout2+1]
		fout[Fout2+0] = fout[fout_ptr+0] - t_r
		fout[Fout2+1] = fout[fout_ptr+1] - t_i
		fout[fout_ptr+0] += t_r
		fout[fout_ptr+1] += t_i

		// 第二组蝶形
		t_r = sMul(fout[Fout2+2]+fout[Fout2+3], tw)
		t_i = sMul(fout[Fout2+3]-fout[Fout2+2], tw)
		fout[Fout2+2] = fout[fout_ptr+2] - t_r
		fout[Fout2+3] = fout[fout_ptr+3] - t_i
		fout[fout_ptr+2] += t_r
		fout[fout_ptr+3] += t_i

		// 第三组蝶形
		t_r = fout[Fout2+5]
		t_i = -fout[Fout2+4]
		fout[Fout2+4] = fout[fout_ptr+4] - t_r
		fout[Fout2+5] = fout[fout_ptr+5] - t_i
		fout[fout_ptr+4] += t_r
		fout[fout_ptr+5] += t_i

		// 第四组蝶形
		t_r = sMul(fout[Fout2+7]-fout[Fout2+6], tw)
		t_i = sMul(-fout[Fout2+7]-fout[Fout2+6], tw)
		fout[Fout2+6] = fout[fout_ptr+6] - t_r
		fout[Fout2+7] = fout[fout_ptr+7] - t_i
		fout[fout_ptr+6] += t_r
		fout[fout_ptr+7] += t_i
	}
}

// kf_bfly4 基4蝴蝶运算
func kf_bfly4(fout []int, fout_ptr, fstride int, st *FFTState, m, N, mm int) {
	if m == 1 {
		// 简化情况: 所有旋转因子为1
		for i := 0; i < N; i++ {
			fp := fout_ptr + i*8
			scratch0 := fout[fp+0] - fout[fp+4]
			scratch1 := fout[fp+1] - fout[fp+5]
			fout[fp+0] += fout[fp+4]
			fout[fp+1] += fout[fp+5]
			scratch2 := fout[fp+2] + fout[fp+6]
			scratch3 := fout[fp+3] + fout[fp+7]
			fout[fp+4] = fout[fp+0] - scratch2
			fout[fp+5] = fout[fp+1] - scratch3
			fout[fp+0] += scratch2
			fout[fp+1] += scratch3
			scratch2 = fout[fp+2] - fout[fp+6]
			scratch3 = fout[fp+3] - fout[fp+7]
			fout[fp+2] = scratch0 + scratch3
			fout[fp+3] = scratch1 - scratch2
			fout[fp+6] = scratch0 - scratch3
			fout[fp+7] = scratch1 + scratch2
		}
	} else {
		// 完整基4蝴蝶
		Fout_beg := fout_ptr
		for i := 0; i < N; i++ {
			fout_ptr = Fout_beg + 2*i*mm
			tw1, tw2, tw3 := 0, 0, 0
			m1 := fout_ptr + 2*m
			m2 := fout_ptr + 4*m
			m3 := fout_ptr + 6*m

			for j := 0; j < m; j++ {
				// 加载输入和旋转因子
				f0 := fout[fout_ptr : fout_ptr+2]
				f1 := fout[m1 : m1+2]
				f2 := fout[m2 : m2+2]
				f3 := fout[m3 : m3+2]
				t1 := st.Twiddles[tw1 : tw1+2]
				t2 := st.Twiddles[tw2 : tw2+2]
				t3 := st.Twiddles[tw3 : tw3+2]

				// 计算蝶形
				scratch0 := sMul(f1[0], t1[0]) - sMul(f1[1], t1[1])
				scratch1 := sMul(f1[0], t1[1]) + sMul(f1[1], t1[0])
				scratch2 := sMul(f2[0], t2[0]) - sMul(f2[1], t2[1])
				scratch3 := sMul(f2[0], t2[1]) + sMul(f2[1], t2[0])
				scratch4 := sMul(f3[0], t3[0]) - sMul(f3[1], t3[1])
				scratch5 := sMul(f3[0], t3[1]) + sMul(f3[1], t3[0])

				// 组合中间结果
				scratch10 := f0[0] - scratch2
				scratch11 := f0[1] - scratch3
				f0[0] += scratch2
				f0[1] += scratch3

				scratch6 := scratch0 + scratch4
				scratch7 := scratch1 + scratch5
				scratch8 := scratch0 - scratch4
				scratch9 := scratch1 - scratch5

				// 更新输出
				f2[0] = f0[0] - scratch6
				f2[1] = f0[1] - scratch7
				f0[0] += scratch6
				f0[1] += scratch7
				f1[0] = scratch10 + scratch9
				f1[1] = scratch11 - scratch8
				f3[0] = scratch10 - scratch9
				f3[1] = scratch11 + scratch8

				// 更新指针和索引
				fout_ptr += 2
				m1 += 2
				m2 += 2
				m3 += 2
				tw1 += 2 * fstride
				tw2 += 4 * fstride
				tw3 += 6 * fstride
			}
		}
	}
}

// kf_bfly3 基3蝴蝶运算
func kf_bfly3(fout []int, fout_ptr, fstride int, st *FFTState, m, N, mm int) {
	Fout_beg := fout_ptr
	m1, m2 := 2*m, 4*m
	ya_tw := Q15_28378 // 预计算的旋转因子

	for i := 0; i < N; i++ {
		fout_ptr = Fout_beg + 2*i*mm
		tw1, tw2 := 0, 0

		for k := 0; k < m; k++ {
			// 加载输入和旋转因子
			f0 := fout[fout_ptr : fout_ptr+2]
			f1 := fout[fout_ptr+m1 : fout_ptr+m1+2]
			f2 := fout[fout_ptr+m2 : fout_ptr+m2+2]
			t1 := st.Twiddles[tw1 : tw1+2]

			// 计算蝶形
			scratch2 := sMul(f1[0], t1[0]) - sMul(f1[1], t1[1])
			scratch3 := sMul(f1[0], t1[1]) + sMul(f1[1], t1[0])
			scratch4 := sMul(f2[0], t1[tw2*2]) - sMul(f2[1], t1[tw2*2+1]) // 简化访问
			scratch5 := sMul(f2[0], t1[tw2*2+1]) + sMul(f2[1], t1[tw2*2]) // 简化访问

			scratch6 := scratch2 + scratch4
			scratch7 := scratch3 + scratch5
			scratch0 := scratch2 - scratch4
			scratch1 := scratch3 - scratch5

			// 更新输出
			f1[0] = f0[0] - halfOf(scratch6)
			f1[1] = f0[1] - halfOf(scratch7)
			f0[0] += scratch6
			f0[1] += scratch7

			scratch0 = sMul(scratch0, ya_tw)
			scratch1 = sMul(scratch1, ya_tw)

			f2[0] = f1[0] + scratch1
			f2[1] = f1[1] - scratch0
			f1[0] -= scratch1
			f1[1] += scratch0

			// 更新指针
			fout_ptr += 2
			tw1 += 2 * fstride
			tw2 += fstride
		}
	}
}

// kf_bfly5 基5蝴蝶运算
func kf_bfly5(fout []int, fout_ptr, fstride int, st *FFTState, m, N, mm int) {
	Fout_beg := fout_ptr
	ya_r, ya_i := 10126, -31164
	yb_r, yb_i := -26510, -19261

	for i := 0; i < N; i++ {
		fout_ptr = Fout_beg + 2*i*mm
		tw1, tw2, tw3, tw4 := 0, 0, 0, 0

		for u := 0; u < m; u++ {
			// 定义指针
			f0 := fout_ptr
			f1 := fout_ptr + 2*m
			f2 := fout_ptr + 4*m
			f3 := fout_ptr + 6*m
			f4 := fout_ptr + 8*m

			// 加载输入并应用旋转因子
			scratch0 := fout[f0]
			scratch1 := fout[f0+1]

			scratch2 := sMul(fout[f1], st.Twiddles[tw1]) - sMul(fout[f1+1], st.Twiddles[tw1+1])
			scratch3 := sMul(fout[f1], st.Twiddles[tw1+1]) + sMul(fout[f1+1], st.Twiddles[tw1])

			scratch4 := sMul(fout[f2], st.Twiddles[tw2]) - sMul(fout[f2+1], st.Twiddles[tw2+1])
			scratch5 := sMul(fout[f2], st.Twiddles[tw2+1]) + sMul(fout[f2+1], st.Twiddles[tw2])

			scratch6 := sMul(fout[f3], st.Twiddles[tw3]) - sMul(fout[f3+1], st.Twiddles[tw3+1])
			scratch7 := sMul(fout[f3], st.Twiddles[tw3+1]) + sMul(fout[f3+1], st.Twiddles[tw3])

			scratch8 := sMul(fout[f4], st.Twiddles[tw4]) - sMul(fout[f4+1], st.Twiddles[tw4+1])
			scratch9 := sMul(fout[f4], st.Twiddles[tw4+1]) + sMul(fout[f4+1], st.Twiddles[tw4])

			// 组合中间结果
			scratch14 := scratch2 + scratch8
			scratch15 := scratch3 + scratch9
			scratch16 := scratch4 + scratch6
			scratch17 := scratch5 + scratch7
			scratch18 := scratch4 - scratch6
			scratch19 := scratch5 - scratch7
			scratch20 := scratch2 - scratch8
			scratch21 := scratch3 - scratch9

			// 更新主输出
			fout[f0] += scratch14 + scratch16
			fout[f0+1] += scratch15 + scratch17

			// 计算中间值
			scratch10 := scratch0 + sMul(scratch14, ya_r) + sMul(scratch16, yb_r)
			scratch11 := scratch1 + sMul(scratch15, ya_r) + sMul(scratch17, yb_r)
			scratch12 := sMul(scratch21, ya_i) + sMul(scratch19, yb_i)
			scratch13 := -sMul(scratch20, ya_i) - sMul(scratch18, yb_i)

			// 更新其他输出
			fout[f1] = scratch10 - scratch12
			fout[f1+1] = scratch11 - scratch13
			fout[f4] = scratch10 + scratch12
			fout[f4+1] = scratch11 + scratch13

			scratch22 := scratch0 + sMul(scratch14, yb_r) + sMul(scratch16, ya_r)
			scratch23 := scratch1 + sMul(scratch15, yb_r) + sMul(scratch17, ya_r)
			scratch24 := -sMul(scratch21, yb_i) + sMul(scratch19, ya_i)
			scratch25 := sMul(scratch20, yb_i) - sMul(scratch18, ya_i)

			fout[f2] = scratch22 + scratch24
			fout[f2+1] = scratch23 + scratch25
			fout[f3] = scratch22 - scratch24
			fout[f3+1] = scratch23 - scratch25

			// 更新位置指针
			fout_ptr += 2
			tw1 += 2 * fstride
			tw2 += 4 * fstride
			tw3 += 6 * fstride
			tw4 += 8 * fstride
		}
	}
}

// OpusFFTImpl FFT核心实现
func OpusFFTImpl(st *FFTState, fout []int, fout_ptr int) {
	fstride := make([]int, MAXFACTORS+1)
	fstride[0] = 1
	L := 0
	m := 1

	// 计算步长因子
	for m != 1 {
		p := st.Factors[2*L]
		m = st.Factors[2*L+1]
		fstride[L+1] = fstride[L] * p
		L++
	}

	// 确定shift值
	shift := st.Shift
	if shift < 0 {
		shift = 0
	}

	// 反向遍历因子
	m = st.Factors[2*L-1]
	for i := L - 1; i >= 0; i-- {
		var m2 int
		if i != 0 {
			m2 = st.Factors[2*i-1]
		} else {
			m2 = 1
		}

		// 根据基数选择蝴蝶操作
		switch p := st.Factors[2*i]; p {
		case 2:
			kf_bfly2(fout, fout_ptr, m, fstride[i])
		case 4:
			kf_bfly4(fout, fout_ptr, fstride[i]<<shift, st, m, fstride[i], m2)
		case 3:
			kf_bfly3(fout, fout_ptr, fstride[i]<<shift, st, m, fstride[i], m2)
		case 5:
			kf_bfly5(fout, fout_ptr, fstride[i]<<shift, st, m, fstride[i], m2)
		}
		m = m2
	}
}

// OpusFFT 执行正向FFT变换
func OpusFFT(st *FFTState, fin, fout []int) {
	scale := st.Scale
	scale_shift := st.ScaleShift - 1

	// 应用缩放和位反转
	for i := 0; i < st.Nfft; i++ {
		rev := st.Bitrev[i]
		fout[2*rev] = SHR32(MULT16_32_Q15(scale, fin[2*i]), scale_shift)
		fout[2*rev+1] = SHR32(MULT16_32_Q15(scale, fin[2*i+1]), scale_shift)
	}

	// 执行FFT核心计算
	OpusFFTImpl(st, fout, 0)
}
