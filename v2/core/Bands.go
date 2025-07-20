package opus


type Bands struct{}

func (b *Bands) Hysteresis_decision(
	val int,
	thresholds []int,
	hysteresis []int,
	N int,
	prev int,
) int {
	i := 0
	for ; i < N; i++ {
		if val < thresholds[i] {
			break
		}
	}

	if i > prev && val < thresholds[prev]+hysteresis[prev] {
		i = prev
	}

	if i < prev && val > thresholds[prev-1]-hysteresis[prev-1] {
		i = prev
	}

	return i
}

func (b *Bands) Celt_lcg_rand(seed int) int {
	return 1664525*seed + 1013904223
}

func (b *Bands) Bitexact_cos(x int) int {
	tmp := 0
	x2 := 0
	tmp = (4096 + x*x) >> 13
	Inlines.OpusAssert(tmp <= 32767)
	x2 = tmp
	x2 = (32767 - x2) + Inlines.FRAC_MUL16(x2, (-7651 + Inlines.FRAC_MUL16(x2, (8277 + Inlines.FRAC_MUL16(-626, x2))))
	Inlines.OpusAssert(x2 <= 32766)
	return 1 + x2
}

func (b *Bands) Bitexact_log2tan(isin, icos int) int {
	lc := Inlines.EC_ILOG(int64(icos))
	ls := Inlines.EC_ILOG(int64(isin))
	icos <<= 15 - lc
	isin <<= 15 - ls
	return (ls-lc)*(1<<11) +
		Inlines.FRAC_MUL16(isin, Inlines.FRAC_MUL16(isin, -2597)+7932) -
		Inlines.FRAC_MUL16(icos, Inlines.FRAC_MUL16(icos, -2597)+7932)
}

func (b *Bands) Compute_band_energies(
	m *CeltMode,
	X [][]int,
	bandE [][]int,
	end int,
	C int,
	LM int,
) {
	i, c, N := 0, 0, 0
	eBands := m.eBands
	N = m.shortMdctSize << LM
	c = 0

	for c < C {
		for i = 0; i < end; i++ {
			j := 0
			maxval := 0
			sum := 0
			maxval = Inlines.Celt_maxabs32(X[c], eBands[i]<<LM, (eBands[i+1]-eBands[i])<<LM)
			if maxval > 0 {
				shift := Inlines.Celt_ilog2(maxval) - 14 + ((m.logN[i]>>EntropyCoder.BITRES + LM + 1) >> 1)
				j = eBands[i] << LM
				if shift > 0 {
					for j < eBands[i+1]<<LM {
						sum = Inlines.MAC16_16(sum, Inlines.EXTRACT16(Inlines.SHR32(X[c][j], shift)), Inlines.EXTRACT16(Inlines.SHR32(X[c][j], shift)))
						j++
					}
				} else {
					for j < eBands[i+1]<<LM {
						sum = Inlines.MAC16_16(sum, Inlines.EXTRACT16(Inlines.SHL32(X[c][j], -shift)), Inlines.EXTRACT16(Inlines.SHL32(X[c][j], -shift)))
						j++
					}
				}
				bandE[c][i] = CeltConstants.EPSILON + Inlines.VSHR32(Inlines.Celt_sqrt(sum), -shift)
			} else {
				bandE[c][i] = CeltConstants.EPSILON
			}
		}
		c++
	}
}

func (b *Bands) Normalise_bands(
	m *CeltMode,
	freq [][]int,
	X [][]int,
	bandE [][]int,
	end int,
	C int,
	M int,
) {
	i, c := 0, 0
	eBands := m.eBands
	c = 0

	for c < C {
		i = 0
		for i < end {
			g := 0
			j, shift := 0, 0
			E := 0
			shift = Inlines.Celt_zlog2(bandE[c][i]) - 13
			E = Inlines.VSHR32(bandE[c][i], shift)
			g = Inlines.EXTRACT16(Inlines.Celt_rcp(Inlines.SHL32(E, 3)))
			j = M * eBands[i]
			for j < M*eBands[i+1] {
				X[c][j] = Inlines.MULT16_16_Q15(Inlines.VSHR32(freq[c][j], shift-1), g)
				j++
			}
			i++
		}
		c++
	}
}

func (b *Bands) Denormalise_bands(
	m *CeltMode,
	X []int,
	freq []int,
	freq_ptr int,
	bandLogE []int,
	bandLogE_ptr int,
	start int,
	end int,
	M int,
	downsample int,
	silence int,
) {
	i, N := 0, 0
	bound := 0
	f := 0
	x := 0
	eBands := m.eBands
	N = M * m.shortMdctSize
	bound = M * eBands[end]
	if downsample != 1 {
		bound = Inlines.IMIN(bound, N/downsample)
	}
	if silence != 0 {
		bound = 0
		start = 0
		end = 0
	}
	f = freq_ptr
	x = M * eBands[start]

	for i = 0; i < M*eBands[start]; i++ {
		freq[f] = 0
		f++
	}

	for i = start; i < end; i++ {
		j, band_end := 0, 0
		g := 0
		lg := 0
		shift := 0

		j = M * eBands[i]
		band_end = M * eBands[i+1]
		lg = Inlines.ADD16(bandLogE[bandLogE_ptr+i], Inlines.SHL16(CeltTables.EMeans[i], 6))

		shift = 16 - (lg >> CeltConstants.DB_SHIFT)
		if shift > 31 {
			shift = 0
			g = 0
		} else {
			g = Inlines.Celt_exp2_frac(lg & ((1 << CeltConstants.DB_SHIFT) - 1))
		}
		if shift < 0 {
			if shift < -2 {
				g = 32767
				shift = -2
			}
			for j < band_end {
				freq[f] = Inlines.SHR32(Inlines.MULT16_16(X[x], g), -shift)
				j++
				x++
				f++
			}
		} else {
			for j < band_end {
				freq[f] = Inlines.SHR32(Inlines.MULT16_16(X[x], g), shift)
				j++
				x++
				f++
			}
		}
	}

	Inlines.OpusAssert(start <= end)
	Arrays.MemSetWithOffset(freq, 0, freq_ptr+bound, N-bound)
}

func (b *Bands) Anti_collapse(
	m *CeltMode,
	X_ [][]int,
	collapse_masks []int16,
	LM int,
	C int,
	size int,
	start int,
	end int,
	logE []int,
	prev1logE []int,
	prev2logE []int,
	pulses []int,
	seed int,
) {
	c, i, j, k := 0, 0, 0, 0
	for i = start; i < end; i++ {
		N0 := 0
		thresh, sqrt_1 := 0, 0
		depth := 0
		shift := 0
		thresh32 := 0

		N0 = m.eBands[i+1] - m.eBands[i]
		Inlines.OpusAssert(pulses[i] >= 0)
		depth = Inlines.Celt_udiv(1+pulses[i], m.eBands[i+1]-m.eBands[i]) >> LM

		thresh32 = Inlines.SHR32(Inlines.Celt_exp2(0-Inlines.SHL16(depth, 10-EntropyCoder.BITRES)), 1)
		thresh = Inlines.MULT16_32_Q15(int16(0.5 * 32767+0.5), Inlines.MIN32(32767, thresh32))
		{
			t := N0 << LM
			shift = Inlines.Celt_ilog2(t) >> 1
			t = Inlines.SHL32(t, (7-shift)<<1)
			sqrt_1 = Inlines.Celt_rsqrt_norm(t)
		}

		c = 0
		for c < C {
			X := 0
			prev1 := 0
			prev2 := 0
			Ediff := 0
			r := 0
			renormalize := 0
			prev1 = prev1logE[c*m.nbEBands+i]
			prev2 = prev2logE[c*m.nbEBands+i]
			if C == 1 {
				prev1 = Inlines.MAX16(prev1, prev1logE[m.nbEBands+i])
				prev2 = Inlines.MAX16(prev2, prev2logE[m.nbEBands+i])
			}
			Ediff = Inlines.EXTEND32(logE[c*m.nbEBands+i]) - Inlines.EXTEND32(Inlines.MIN16(prev1, prev2))
			Ediff = Inlines.MAX32(0, Ediff)

			if Ediff < 16384 {
				r32 := Inlines.SHR32(Inlines.Celt_exp2(int16(0-Inlines.EXTRACT16(Ediff))), 1)
				r = 2 * Inlines.MIN16(16383, r32)
			} else {
				r = 0
			}
			if LM == 3 {
				r = Inlines.MULT16_16_Q14(23170, Inlines.MIN32(23169, r))
			}
			r = Inlines.SHR16(Inlines.MIN16(thresh, r), 1)
			r = Inlines.SHR32(Inlines.MULT16_16_Q15(sqrt_1, r), shift)

			X = m.eBands[i] << LM
			for k = 0; k < 1<<LM; k++ {
				if (collapse_masks[i*C+c] & (1 << k)) == 0 {
					Xk := X + k
					for j = 0; j < N0; j++ {
						seed = b.Celt_lcg_rand(seed)
						if (seed & 0x8000) != 0 {
							X_[c][Xk+j<<LM] = r
						} else {
							X_[c][Xk+j<<LM] = -r
						}
					}
					renormalize = 1
				}
			}
			if renormalize != 0 {
				VQ.Renormalise_vector(X_[c], X, N0<<LM, CeltConstants.Q15ONE)
			}
			c++
		}
	}
}

func (b *Bands) Intensity_stereo(
	m *CeltMode,
	X []int,
	X_ptr int,
	Y []int,
	Y_ptr int,
	bandE [][]int,
	bandID int,
	N int,
) {
	i := bandID
	j := 0
	a1, a2 := 0, 0
	left, right := 0, 0
	norm := 0
	shift := Inlines.Celt_zlog2(Inlines.MAX32(bandE[0][i], bandE[1][i])) - 13
	left = Inlines.VSHR32(bandE[0][i], shift)
	right = Inlines.VSHR32(bandE[1][i], shift)
	norm = CeltConstants.EPSILON + Inlines.Celt_sqrt(CeltConstants.EPSILON+Inlines.MULT16_16(left, left)+Inlines.MULT16_16(right, right))
	a1 = Inlines.DIV32_16(Inlines.SHL32(left, 14), norm)
	a2 = Inlines.DIV32_16(Inlines.SHL32(right, 14), norm)
	for j = 0; j < N; j++ {
		r, l := Y[Y_ptr+j], X[X_ptr+j]
		X[X_ptr+j] = Inlines.EXTRACT16(Inlines.SHR32(Inlines.MAC16_16(Inlines.MULT16_16(a1, l), a2, r), 14))
	}
}

func (b *Bands) Stereo_split(
	X []int,
	X_ptr int,
	Y []int,
	Y_ptr int,
	N int,
) {
	j := 0
	for j = 0; j < N; j++ {
		r, l := Y[Y_ptr+j], X[X_ptr+j]
		l = Inlines.MULT16_16(int16(0.70710678 * 32767+0.5), l)
		r = Inlines.MULT16_16(int16(0.70710678 * 32767+0.5), r)
		X[X_ptr+j] = Inlines.EXTRACT16(Inlines.SHR32(Inlines.ADD32(l, r), 15))
		Y[Y_ptr+j] = Inlines.EXTRACT16(Inlines.SHR32(Inlines.SUB32(r, l), 15))
	}
}

func (b *Bands) Stereo_merge(
	X []int,
	X_ptr int,
	Y []int,
	Y_ptr int,
	mid int,
	N int,
) {
	j := 0
	xp := &BoxedValueInt{Val: 0}
	side := &BoxedValueInt{Val: 0}
	El, Er := 0, 0
	mid2 := 0
	kl, kr := 0, 0
	t, lgain, rgain := 0, 0, 0

	Kernels.Dual_inner_prod(Y, Y_ptr, X, X_ptr, Y, Y_ptr, N, xp, side)
	xp.Val = Inlines.MULT16_32_Q15(int16(mid), xp.Val)
	mid2 = Inlines.SHR16(int16(mid), 1)
	El = Inlines.MULT16_16(mid2, mid2) + side.Val - 2*xp.Val
	Er = Inlines.MULT16_16(mid2, mid2) + side.Val + 2*xp.Val
	if Er < int(0.5+6e-4*(1<<28)) || El < int(0.5+6e-4*(1<<28)) {
		copy(Y[Y_ptr:Y_ptr+N], X[X_ptr:X_ptr+N])
		return
	}

	kl = Inlines.Celt_ilog2(El) >> 1
	kr = Inlines.Celt_ilog2(Er) >> 1
	t = Inlines.VSHR32(El, (kl-7)<<1)
	lgain = Inlines.Celt_rsqrt_norm(t)
	t = Inlines.VSHR32(Er, (kr-7)<<1)
	rgain = Inlines.Celt_rsqrt_norm(t)

	if kl < 7 {
		kl = 7
	}
	if kr < 7 {
		kr = 7
	}

	for j = 0; j < N; j++ {
		r, l := Y[Y_ptr+j], Inlines.MULT16_16_P15(int16(mid), X[X_ptr+j])
		X[X_ptr+j] = Inlines.EXTRACT16(Inlines.PSHR32(Inlines.MULT16_16(int16(lgain), Inlines.SUB16(l, r)), kl+1))
		Y[Y_ptr+j] = Inlines.EXTRACT16(Inlines.PSHR32(Inlines.MULT16_16(int16(rgain), Inlines.ADD16(l, r)), kr+1))
	}
}

func (b *Bands) Spreading_decision(
	m *CeltMode,
	X [][]int,
	average *BoxedValueInt,
	last_decision int,
	hf_average *BoxedValueInt,
	tapset_decision *BoxedValueInt,
	update_hf int,
	end int,
	C int,
	M int,
) int {
	i, c := 0, 0
	sum, nbBands := 0, 0
	eBands := m.eBands
	decision := 0
	hf_sum := 0

	Inlines.OpusAssert(end > 0)

	if M*(eBands[end]-eBands[end-1]) <= 8 {
		return Spread.SPREAD_NONE
	}

	c = 0
	for c < C {
		for i = 0; i < end; i++ {
			j, N, tmp := 0, 0, 0
			tcount := [3]int{0, 0, 0}
			x := X[c]
			x_ptr := M * eBands[i]
			N = M * (eBands[i+1] - eBands[i])
			if N <= 8 {
				continue
			}
			for j = x_ptr; j < x_ptr+N; j++ {
				x2N := Inlines.MULT16_16(Inlines.MULT16_16_Q15(x[j], x[j]), N)
				if x2N < int(0.25 * 8192+0.5) {
					tcount[0]++
				}
				if x2N < int(0.0625 * 8192+0.5) {
					tcount[1]++
				}
				if x2N < int(0.015625 * 8192+0.5) {
					tcount[2]++
				}
			}
			if i > m.nbEBands-4 {
				hf_sum += Inlines.Celt_udiv(32*(tcount[1]+tcount[0]), N)
			}
			tmp = 0
			if 2*tcount[2] >= N {
				tmp++
			}
			if 2*tcount[1] >= N {
				tmp++
			}
			if 2*tcount[0] >= N {
				tmp++
			}
			sum += tmp * 256
			nbBands++
		}
		c++
	}

	if update_hf != 0 {
		if hf_sum != 0 {
			hf_sum = Inlines.Celt_udiv(hf_sum, C*(4-m.nbEBands+end))
		}
		hf_average.Val = (hf_average.Val + hf_sum) >> 1
		hf_sum = hf_average.Val
		if tapset_decision.Val == 2 {
			hf_sum += 4
		} else if tapset_decision.Val == 0 {
			hf_sum -= 4
		}
		if hf_sum > 22 {
			tapset_decision.Val = 2
		} else if hf_sum > 18 {
			tapset_decision.Val = 1
		} else {
			tapset_decision.Val = 0
		}
	}

	Inlines.OpusAssert(nbBands > 0)
	sum = Inlines.Celt_udiv(sum, nbBands)
	sum = (sum + average.Val) >> 1
	average.Val = sum
	sum = (3*sum + ((3-last_decision)<<7+64) + 2) >> 2
	if sum < 80 {
		decision = Spread.SPREAD_AGGRESSIVE
	} else if sum < 256 {
		decision = Spread.SPREAD_NORMAL
	} else if sum < 384 {
		decision = Spread.SPREAD_LIGHT
	} else {
		decision = Spread.SPREAD_NONE
	}
	return decision
}

func (b *Bands) Deinterleave_hadamard(
	X []int,
	X_ptr int,
	N0 int,
	stride int,
	hadamard int,
) {
	i, j := 0, 0
	N := N0 * stride
	tmp := make([]int, N)

	Inlines.OpusAssert(stride > 0)
	if hadamard != 0 {
		ordery := stride - 2
		for i = 0; i < stride; i++ {
			for j = 0; j < N0; j++ {
				tmp[CeltTables.Ordery_table[ordery+i]*N0+j] = X[X_ptr+j*stride+i]
			}
		}
	} else {
		for i = 0; i < stride; i++ {
			for j = 0; j < N0; j++ {
				tmp[i*N0+j] = X[X_ptr+j*stride+i]
			}
		}
	}
	copy(X[X_ptr:X_ptr+N], tmp)
}

func (b *Bands) Interleave_hadamard(
	X []int,
	X_ptr int,
	N0 int,
	stride int,
	hadamard int,
) {
	i, j := 0, 0
	N := N0 * stride
	tmp := make([]int, N)

	if hadamard != 0 {
		ordery := stride - 2
		for i = 0; i < stride; i++ {
			for j = 0; j < N0; j++ {
				tmp[j*stride+i] = X[X_ptr+CeltTables.Ordery_table[ordery+i]*N0+j]
			}
		}
	} else {
		for i = 0; i < stride; i++ {
			for j = 0; j < N0; j++ {
				tmp[j*stride+i] = X[X_ptr+i*N0+j]
			}
		}
	}
	copy(X[X_ptr:X_ptr+N], tmp)
}

func (b *Bands) Haar1(
	X []int,
	X_ptr int,
	N0 int,
	stride int,
) {
	i, j := 0, 0
	N0 >>= 1
	for i = 0; i < stride; i++ {
		for j = 0; j < N0; j++ {
			tmpidx := X_ptr + i + stride*2*j
			tmp1 := Inlines.MULT16_16(int16(0.70710678 * 32767+0.5), X[tmpidx])
			tmp2 := Inlines.MULT16_16(int16(0.70710678 * 32767+0.5), X[tmpidx+stride])
			X[tmpidx] = Inlines.EXTRACT16(Inlines.PSHR32(Inlines.ADD32(tmp1, tmp2), 15))
			X[tmpidx+stride] = Inlines.EXTRACT16(Inlines.PSHR32(Inlines.SUB32(tmp1, tmp2), 15))
		}
	}
}

func (b *Bands) Haar1ZeroOffset(
	X []int,
	N0 int,
	stride int,
) {
	i, j := 0, 0
	N0 >>= 1
	for i = 0; i < stride; i++ {
		for j = 0; j < N0; j++ {
			tmpidx := i + stride*2*j
			tmp1 := Inlines.MULT16_16(int16(0.70710678 * 32767+0.5), X[tmpidx])
			tmp2 := Inlines.MULT16_16(int16(0.70710678 * 32767+0.5), X[tmpidx+stride])
			X[tmpidx] = Inlines.EXTRACT16(Inlines.PSHR32(Inlines.ADD32(tmp1, tmp2), 15))
			X[tmpidx+stride] = Inlines.EXTRACT16(Inlines.PSHR32(Inlines.SUB32(tmp1, tmp2), 15))
		}
	}
}

func (b *Bands) Compute_qn(
	N int,
	bits int,
	offset int,
	pulse_cap int,
	stereo int,
) int {
	exp2_table8 := []int16{16384, 17866, 19483, 21247, 23170, 25267, 27554, 30048}
	qn, qb := 0, 0
	N2 := 2*N - 1
	if stereo != 0 && N == 2 {
		N2--
	}
	qb = Inlines.Celt_sudiv(bits+N2*offset, N2)
	qb = Inlines.IMIN(bits-pulse_cap-(4<<EntropyCoder.BITRES), qb)
	qb = Inlines.IMIN(8<<EntropyCoder.BITRES, qb)
	if qb < (1<<EntropyCoder.BITRES)>>1 {
		qn = 1
	} else {
		qn = int(exp2_table8[qb&0x7]) >> (14 - (qb >> EntropyCoder.BITRES))
		qn = (qn+1)>>1<<1
	}
	Inlines.OpusAssert(qn <= 256)
	return qn
}

type Band_ctx struct {
	Encode        int
	M             *CeltMode
	I             int
	Intensity     int
	Spread        int
	Tf_change     int
	Ec            *EntropyCoder.EC
	Remaining_bits int
	BandE         [][]int
	Seed          int
}

type Split_ctx struct {
	Inv    int
	Imid   int
	Iside  int
	Delta  int
	Itheta int
	Qalloc int
}

func (b *Bands) Compute_theta(
	ctx *Band_ctx,
	sctx *Split_ctx,
	X []int,
	X_ptr int,
	Y []int,
	Y_ptr int,
	N int,
	bits *BoxedValueInt,
	B int,
	B0 int,
	LM int,
	stereo int,
	fill *BoxedValueInt,
) {
	qn := 0
	itheta := 0
	delta := 0
	imid, iside := 0, 0
	qalloc := 0
	pulse_cap := 0
	offset := 0
	tell := 0
	inv := 0
	encode := ctx.Encode
	m := ctx.M
	i := ctx.I
	intensity := ctx.Intensity
	ec := ctx.Ec
	bandE := ctx.BandE

	pulse_cap = m.LogN[i] + LM*(1<<EntropyCoder.BITRES)
	offset = (pulse_cap >> 1) - Ternary(stereo != 0 && N == 2, CeltConstants.QTHETA_OFFSET_TWOPHASE, CeltConstants.QTHETA_OFFSET)
	qn = b.Compute_qn(N, bits.Val, offset, pulse_cap, stereo)
	if stereo != 0 && i >= intensity {
		qn = 1
	}

	if encode != 0 {
		itheta = VQ.Stereo_itheta(X, X_ptr, Y, Y_ptr, stereo, N)
	}

	tell = int(ec.Tell_frac())

	if qn != 1 {
		if encode != 0 {
			itheta = (itheta*qn + 8192) >> 14
		}

		if stereo != 0 && N > 2 {
			p0 := 3
			x := itheta
			x0 := qn / 2
			ft := Inlines.CapToUInt32(p0*(x0+1) + x0)
			if encode != 0 {
				if x <= x0 {
					ec.Encode(p0*x, p0*(x+1), ft)
				} else {
					ec.Encode((x-1-x0)+(x0+1)*p0, (x-x0)+(x0+1)*p0, ft)
				}
			} else {
				fs := int(ec.Decode(ft))
				if fs < (x0+1)*p0 {
					x = fs / p0
				} else {
					x = x0 + 1 + (fs - (x0+1)*p0)
				}
				if x <= x0 {
					ec.Dec_update(p0*x, p0*(x+1), ft)
				} else {
					ec.Dec_update((x-1-x0)+(x0+1)*p0, (x-x0)+(x0+1)*p0, ft)
				}
				itheta = x
			}
		} else if B0 > 1 || stereo != 0 {
			if encode != 0 {
				ec.Enc_uint(itheta, qn+1)
			} else {
				itheta = int(ec.Dec_uint(qn + 1))
			}
		} else {
			fs := 1
			ft := ((qn >> 1) + 1) * ((qn >> 1) + 1)
			if encode != 0 {
				fl := 0
				if itheta <= qn>>1 {
					fs = itheta + 1
					fl = itheta * (itheta + 1) >> 1
				} else {
					fs = qn + 1 - itheta
					fl = ft - (qn+1-itheta)*(qn+2-itheta)>>1
				}
				ec.Encode(fl, fl+fs, ft)
			} else {
				fl := 0
				fm := int(ec.Decode(ft))
				if fm < ((qn>>1)*((qn>>1)+1))>>1 {
					itheta = (Inlines.Isqrt32(8*fm+1) - 1) >> 1
					fs = itheta + 1
					fl = itheta * (itheta + 1) >> 1
				} else {
					itheta = (2*(qn+1) - Inlines.Isqrt32(8*(ft-fm-1)+1)) >> 1
					fs = qn + 1 - itheta
					fl = ft - (qn+1-itheta)*(qn+2-itheta)>>1
				}
				ec.Dec_update(fl, fl+fs, ft)
			}
		}
		Inlines.OpusAssert(itheta >= 0)
		itheta = Inlines.Celt_udiv(itheta*16384, qn)
		if encode != 0 && stereo != 0 {
			if itheta == 0 {
				b.Intensity_stereo(m, X, X_ptr, Y, Y_ptr, bandE, i, N)
			} else {
				b.Stereo_split(X, X_ptr, Y, Y_ptr, N)
			}
		}
	} else if stereo != 0 {
		if encode != 0 {
			inv = Ternary(itheta > 8192, 1, 0)
			if inv != 0 {
				for j := 0; j < N; j++ {
					Y[Y_ptr+j] = -Y[Y_ptr+j]
				}
			}
			b.Intensity_stereo(m, X, X_ptr, Y, Y_ptr, bandE, i, N)
		}
		if bits.Val > 2<<EntropyCoder.BITRES && ctx.Remaining_bits > 2<<EntropyCoder.BITRES {
			if encode != 0 {
				ec.Enc_bit_logp(inv, 2)
			} else {
				inv = ec.Dec_bit_logp(2)
			}
		} else {
			inv = 0
		}
		itheta = 0
	}
	qalloc = int(ec.Tell_frac()) - tell
	bits.Val -= qalloc

	if itheta == 0 {
		imid = 32767
		iside = 0
		fill.Val &= (1 << B) - 1
		delta = -16384
	} else if itheta == 16384 {
		imid = 0
		iside = 32767
		fill.Val &= ((1 << B) - 1) << B
		delta = 16384
	} else {
		imid = b.Bitexact_cos(itheta)
		iside = b.Bitexact_cos(16384 - itheta)
		delta = Inlines.FRAC_MUL16((N-1)<<7, b.Bitexact_log2tan(iside, imid))
	}

	sctx.Inv = inv
	sctx.Imid = imid
	sctx.Iside = iside
	sctx.Delta = delta
	sctx.Itheta = itheta
	sctx.Qalloc = qalloc
}

func (b *Bands) Quant_band_n1(
	ctx *Band_ctx,
	X []int,
	X_ptr int,
	Y []int,
	Y_ptr int,
	bits int,
	lowband_out []int,
	lowband_out_ptr int,
) int {
	resynth := Ternary(ctx.Encode == 0, 1, 0)
	c := 0
	x := X
	x_ptr := X_ptr
	encode := ctx.Encode
	ec := ctx.Ec

	stereo := Ternary(Y != nil, 1, 0)
	for c < 1+stereo {
		sign := 0
		if ctx.Remaining_bits >= 1<<EntropyCoder.BITRES {
			if encode != 0 {
				sign = Ternary(x[x_ptr] < 0, 1, 0)
				ec.Enc_bits(sign, 1)
			} else {
				sign = ec.Dec_bits(1)
			}
			ctx.Remaining_bits -= 1 << EntropyCoder.BITRES
			bits -= 1 << EntropyCoder.BITRES
		}
		if resynth != 0 {
			x[x_ptr] = Ternary(sign != 0, -CeltConstants.NORM_SCALING, CeltConstants.NORM_SCALING)
		}
		x = Y
		x_ptr = Y_ptr
		c++
	}
	if lowband_out != nil {
		lowband_out[lowband_out_ptr] = Inlines.SHR16(X[X_ptr], 4)
	}
	return 1
}

func (b *Bands) Quant_partition(
	ctx *Band_ctx,
	X []int,
	X_ptr int,
	N int,
	bits int,
	B int,
	lowband []int,
	lowband_ptr int,
	LM int,
	gain int,
	fill int,
) int {
	cache_ptr := 0
	q := 0
	curr_bits := 0
	imid, iside := 0, 0
	B0 := B
	mid, side := 0, 0
	cm := 0
	resynth := Ternary(ctx.Encode == 0, 1, 0)
	Y := 0
	encode := ctx.Encode
	m := ctx.M
	i := ctx.I
	spread := ctx.Spread
	ec := ctx.Ec

	cache := m.Cache.Bits
	cache_ptr = m.Cache.Index[(LM+1)*m.nbEBands+i]
	if LM != -1 && bits > cache[cache_ptr]+cache[cache_ptr+1]+12 && N > 2 {
		mbits, sbits, delta := 0, 0, 0
		itheta := 0
		qalloc := 0
		sctx := &Split_ctx{}
		next_lowband2 := 0
		rebalance := 0

		N >>= 1
		Y = X_ptr + N
		LM--
		if B == 1 {
			fill = (fill & 1) | (fill << 1)
		}
		B = (B + 1) >> 1

		boxed_b := &BoxedValueInt{Val: bits}
		boxed_fill := &BoxedValueInt{Val: fill}
		b.Compute_theta(ctx, sctx, X, X_ptr, X, Y, N, boxed_b, B, B0, LM, 0, boxed_fill)
		bits = boxed_b.Val
		fill = boxed_fill.Val

		imid = sctx.Imid
		iside = sctx.Iside
		delta = sctx.Delta
		itheta = sctx.Itheta
		qalloc = sctx.Qalloc
		mid = imid
		side = iside

		if B0 > 1 && (itheta&0x3fff) != 0 {
			if itheta > 8192 {
				delta -= delta >> (4 - LM)
			} else {
				delta = Inlines.IMIN(0, delta+(N<<EntropyCoder.BITRES>>(5-LM)))
			}
		}
		mbits = Inlines.IMAX(0, Inlines.IMIN(bits, (bits-delta)>>1))
		sbits = bits - mbits
		ctx.Remaining_bits -= qalloc

		if lowband != nil {
			next_lowband2 = lowband_ptr + N
		}

		rebalance = ctx.Remaining_bits
		if mbits >= sbits {
			cm = b.Quant_partition(ctx, X, X_ptr, N, mbits, B, lowband, lowband_ptr, LM, Inlines.MULT16_16_P15(int16(gain), int16(mid)), fill)
			rebalance = mbits - (rebalance - ctx.Remaining_bits)
			if rebalance > 3<<EntropyCoder.BITRES && itheta != 0 {
				sbits += rebalance - (3 << EntropyCoder.BITRES)
			}
			cm |= b.Quant_partition(ctx, X, Y, N, sbits, B, lowband, next_lowband2, LM, Inlines.MULT16_16_P15(int16(gain), int16(side)), fill>>B) << (B0 >> 1)
		} else {
			cm = b.Quant_partition(ctx, X, Y, N, sbits, B, lowband, next_lowband2, LM, Inlines.MULT16_16_P15(int16(gain), int16(side)), fill>>B) << (B0 >> 1)
			rebalance = sbits - (rebalance - ctx.Remaining_bits)
			if rebalance > 3<<EntropyCoder.BITRES && itheta != 16384 {
				mbits += rebalance - (3 << EntropyCoder.BITRES)
			}
			cm |= b.Quant_partition(ctx, X, X_ptr, N, mbits, B, lowband, lowband_ptr, LM, Inlines.MULT16_16_P15(int16(gain), int16(mid)), fill)
		}
	} else {
		q = Rate.Bits2pulses(m, i, LM, bits)
		curr_bits = Rate.Pulses2bits(m, i, LM, q)
		ctx.Remaining_bits -= curr_bits

		for ctx.Remaining_bits < 0 && q > 0 {
			ctx.Remaining_bits += curr_bits
			q--
			curr_bits = Rate.Pulses2bits(m, i, LM, q)
			ctx.Remaining_bits -= curr_bits
		}

		if q != 0 {
			K := Rate.Get_pulses(q)
			if encode != 0 {
				cm = VQ.Alg_quant(X, X_ptr, N, K, spread, B, ec)
			} else {
				cm = VQ.Alg_unquant(X, X_ptr, N, K, spread, B, ec, int16(gain))
			}
		} else if resynth != 0 {
			cm_mask := (1 << B) - 1
			fill &= cm_mask
			if fill == 0 {
				Arrays.MemSetWithOffset(X, 0, X_ptr, N)
			} else {
				if lowband == nil {
					for j := 0; j < N; j++ {
						ctx.Seed = b.Celt_lcg_rand(ctx.Seed)
						X[X_ptr+j] = int(ctx.Seed) >> 20
					}
					cm = cm_mask
				} else {
					for j := 0; j < N; j++ {
						tmp := int16(1.0/256 * 1024 + 0.5)
						if (ctx.Seed & 0x8000) != 0 {
							tmp = -tmp
						}
						X[X_ptr+j] = lowband[lowband_ptr+j] + int(tmp)
					}
					cm = fill
				}
				VQ.Renormalise_vector(X, X_ptr, N, int16(gain))
			}
		}
	}
	return cm
}

var bit_interleave_table = []byte{0, 1, 1, 1, 2, 3, 3, 3, 2, 3, 3, 3, 2, 3, 3, 3}
var bit_deinterleave_table = []int{0x00, 0x03, 0x0C, 0x0F, 0x30, 0x33, 0x3C, 0x3F, 0xC0, 0xC3, 0xCC, 0xCF, 0xF0, 0xF3, 0xFC, 0xFF}

func (b *Bands) Quant_band(
	ctx *Band_ctx,
	X []int,
	X_ptr int,
	N int,
	bits int,
	B int,
	lowband []int,
	lowband_ptr int,
	LM int,
	lowband_out []int,
	lowband_out_ptr int,
	gain int,
	lowband_scratch []int,
	lowband_scratch_ptr int,
	fill int,
) int {
	N0 := N
	N_B := N
	N_B0 := 0
	B0 := B
	time_divide := 0
	recombine := 0
	longBlocks := 0
	cm := 0
	resynth := Ternary(ctx.Encode == 0, 1, 0)
	k := 0
	encode := ctx.Encode
	tf_change := ctx.Tf_change

	longBlocks = Ternary(B0 == 1, 1, 0)
	N_B = Inlines.Celt_udiv(N_B, B)

	if N == 1 {
		return b.Quant_band_n1(ctx, X, X_ptr, nil, 0, bits, lowband_out, lowband_out_ptr)
	}

	if tf_change > 0 {
		recombine = tf_change
	}

	if lowband_scratch != nil && lowband != nil && (recombine != 0 || (N_B&1 == 0 && tf_change < 0) || B0 > 1) {
		copy(lowband_scratch[lowband_scratch_ptr:], lowband[lowband_ptr:lowband_ptr+N])
		lowband = lowband_scratch
		lowband_ptr = lowband_scratch_ptr
	}

	for k = 0; k < recombine; k++ {
		if encode != 0 {
			b.Haar1(X, X_ptr, N>>k, 1<<k)
		}
		if lowband != nil {
			b.Haar1(lowband, lowband_ptr, N>>k, 1<<k)
		}
		fill = int(bit_interleave_table[fill&0xF]) | int(bit_interleave_table[fill>>4])<<2
	}
	B >>= recombine
	N_B <<= recombine

	for (N_B&1) == 0 && tf_change < 0 {
		if encode != 0 {
			b.Haar1(X, X_ptr, N_B, B)
		}
		if lowband != nil {
			b.Haar1(lowband, lowband_ptr, N_B, B)
		}
		fill |= fill << B
		B <<= 1
		N_B >>= 1
		time_divide++
		tf_change++
	}
	B0 = B
	N_B0 = N_B

	if B0 > 1 {
		if encode != 0 {
			b.Deinterleave_hadamard(X, X_ptr, N_B>>recombine, B0<<recombine, longBlocks)
		}
		if lowband != nil {
			b.Deinterleave_hadamard(lowband, lowband_ptr, N_B>>recombine, B0<<recombine, longBlocks)
		}
	}

	cm = b.Quant_partition(ctx, X, X_ptr, N, bits, B, lowband, lowband_ptr, LM, gain, fill)

	if resynth != 0 {
		if B0 > 1 {
			b.Interleave_hadamard(X, X_ptr, N_B>>recombine, B0<<recombine, longBlocks)
		}

		N_B = N_B0
		B = B0
		for k = 0; k < time_divide; k++ {
			B >>= 1
			N_B <<= 1
			cm |= cm >> B
			b.Haar1(X, X_ptr, N_B, B)
		}

		for k = 0; k < recombine; k++ {
			cm = bit_deinterleave_table[cm]
			b.Haar1(X, X_ptr, N0>>k, 1<<k)
		}
		B <<= recombine

		if lowband_out != nil {
			n := Inlines.Celt_sqrt(Inlines.SHL32(N0, 22))
			for j := 0; j < N0; j++ {
				lowband_out[lowband_out_ptr+j] = Inlines.MULT16_16_Q15(int16(n), X[X_ptr+j])
			}
		}
		cm &= (1 << B) - 1
	}
	return cm
}

func (b *Bands) Quant_band_stereo(
	ctx *Band_ctx,
	X []int,
	X_ptr int,
	Y []int,
	Y_ptr int,
	N int,
	bits int,
	B int,
	lowband []int,
	lowband_ptr int,
	LM int,
	lowband_out []int,
	lowband_out_ptr int,
	lowband_scratch []int,
	lowband_scratch_ptr int,
	fill int,
) int {
	imid, iside := 0, 0
	inv := 0
	mid, side := 0, 0
	cm := 0
	resynth := Ternary(ctx.Encode == 0, 1, 0)
	mbits, sbits, delta := 0, 0, 0
	itheta := 0
	qalloc := 0
	sctx := &Split_ctx{}
	orig_fill := fill
	encode := ctx.Encode
	ec := ctx.Ec

	if N == 1 {
		return b.Quant_band_n1(ctx, X, X_ptr, Y, Y_ptr, bits, lowband_out, lowband_out_ptr)
	}

	boxed_b := &BoxedValueInt{Val: bits}
	boxed_fill := &BoxedValueInt{Val: fill}
	b.Compute_theta(ctx, sctx, X, X_ptr, Y, Y_ptr, N, boxed_b, B, B, LM, 1, boxed_fill)
	bits = boxed_b.Val
	fill = boxed_fill.Val

	inv = sctx.Inv
	imid = sctx.Imid
	iside = sctx.Iside
	delta = sctx.Delta
	itheta = sctx.Itheta
	qalloc = sctx.Qalloc
	mid = imid
	side = iside

	if N == 2 {
		c := 0
		sign := 0
		var x2, y2 []int
		var x2_ptr, y2_ptr int
		mbits = bits
		sbits = 0
		if itheta != 0 && itheta != 16384 {
			sbits = 1 << EntropyCoder.BITRES
		}
		mbits -= sbits
		c = Ternary(itheta > 8192, 1, 0)
		ctx.Remaining_bits -= qalloc + sbits
		if c != 0 {
			x2 = Y
			x2_ptr = Y_ptr
			y2 = X
			y2_ptr = X_ptr
		} else {
			x2 = X
			x2_ptr = X_ptr
			y2 = Y
			y2_ptr = Y_ptr
		}

		if sbits != 0 {
			if encode != 0 {
				sign = Ternary(X[X_ptr]*Y[Y_ptr+1]-X[X_ptr+1]*Y[Y_ptr] < 0, 1, 0)
				ec.Enc_bits(sign, 1)
			} else {
				sign = ec.Dec_bits(1)
			}
		}
		sign = 1 - 2*sign
		cm = b.Quant_band(ctx, x2, x2_ptr, N, mbits, B, lowband, lowband_ptr,
			LM, lowband_out, lowband_out_ptr, CeltConstants.Q15ONE, lowband_scratch, lowband_scratch_ptr, orig_fill)

		y2[y2_ptr] = -sign * x2[x2_ptr+1]
		y2[y2_ptr+1] = sign * x2[x2_ptr]
		if resynth != 0 {
			tmp := X[X_ptr]
			X[X_ptr] = Inlines.MULT16_16_Q15(int16(mid), X[X_ptr])
			X[X_ptr+1] = Inlines.MULT16_16_Q15(int16(mid), X[X_ptr+1])
			Y[Y_ptr] = Inlines.MULT16_16_Q15(int16(side), Y[Y_ptr])
			Y[Y_ptr+1] = Inlines.MULT16_16_Q15(int16(side), Y[Y_ptr+1])
			X[X_ptr] = Inlines.SUB16(tmp, Y[Y_ptr])
			Y[Y_ptr] = Inlines.ADD16(tmp, Y[Y_ptr])
			tmp = X[X_ptr+1]
			X[X_ptr+1] = Inlines.SUB16(tmp, Y[Y_ptr+1])
			Y[Y_ptr+1] = Inlines.ADD16(tmp, Y[Y_ptr+1])
		}
	} else {
		rebalance := 0
		mbits = Inlines.IMAX(0, Inlines.IMIN(bits, (bits-delta)>>1))
		sbits = bits - mbits
		ctx.Remaining_bits -= qalloc
		rebalance = ctx.Remaining_bits
		if mbits >= sbits {
			cm = b.Quant_band(ctx, X, X_ptr, N, mbits, B,
				lowband, lowband_ptr, LM, lowband_out, lowband_out_ptr,
				CeltConstants.Q15ONE, lowband_scratch, lowband_scratch_ptr, fill)
			rebalance = mbits - (rebalance - ctx.Remaining_bits)
			if rebalance > 3<<EntropyCoder.BITRES && itheta != 0 {
				sbits += rebalance - (3 << EntropyCoder.BITRES)
			}
			cm |= b.Quant_band(ctx, Y, Y_ptr, N, sbits, B,
				nil, 0, LM, nil, 0,
				int16(side), nil, 0, fill>>B)
		} else {
			cm = b.Quant_band(ctx, Y, Y_ptr, N, sbits, B,
				nil, 0, LM, nil, 0,
				int16(side), nil, 0, fill>>B)
			rebalance = sbits - (rebalance - ctx.Remaining_bits)
			if rebalance > 3<<EntropyCoder.BITRES && itheta != 16384 {
				mbits += rebalance - (3 << EntropyCoder.BITRES)
			}
			cm |= b.Quant_band(ctx, X, X_ptr, N, mbits, B,
				lowband, lowband_ptr, LM, lowband_out, lowband_out_ptr,
				CeltConstants.Q15ONE, lowband_scratch, lowband_scratch_ptr, fill)
		}
	}

	if resynth != 0 {
		if N != 2 {
			b.Stereo_merge(X, X_ptr, Y, Y_ptr, mid, N)
		}
		if inv != 0 {
			for j := Y_ptr; j < Y_ptr+N; j++ {
				Y[j] = -Y[j]
			}
		}
	}
	return cm
}

func (b *Bands) Quant_all_bands(
	encode int,
	m *CeltMode,
	start int,
	end int,
	X_ [][]int,
	Y_ [][]int,
	collapse_masks []int16,
	bandE [][]int,
	pulses []int,
	shortBlocks int,
	spread int,
	dual_stereo int,
	intensity int,
	tf_res []int,
	total_bits int,
	balance int,
	ec *EntropyCoder.EC,
	LM int,
	codedBands int,
	seed *BoxedValueInt,
) {
	i := 0
	remaining_bits := 0
	eBands := m.eBands
	norm := make([]int, len(X_))
	norm2 := 0
	lowband_scratch := make([]int, len(X_))
	lowband_scratch_ptr := 0
	B := 0
	M := 0
	lowband_offset := 0
	update_lowband := 1
	C := Ternary(Y_ != nil, 2, 1)
	norm_offset := 0
	resynth := Ternary(encode == 0, 1, 0)
	ctx := &Band_ctx{}

	M = 1 << LM
	B = Ternary(shortBlocks != 0, M, 1)
	norm_offset = M * eBands[start]
	norm2 = M*eBands[m.nbEBands-1] - norm_offset
	lowband_scratch_ptr = M * eBands[m.nbEBands-1]
	lowband_offset = 0
	ctx.BandE = bandE
	ctx.Ec = ec
	ctx.Encode = encode
	ctx.Intensity = intensity
	ctx.M = m
	ctx.Seed = seed.Val
	ctx.Spread = spread
	for i = start; i < end; i++ {
		tell := 0
		bits := 0
		N := 0
		curr_balance := 0
		effective_lowband := -1
		var X, Y []int
		X_ptr, Y_ptr := 0, 0
		last := 0

		ctx.I = i
		last = Ternary(i == end-1, 1, 0)

		X = X_[i]
		X_ptr = M * eBands[i]
		if Y_ != nil {
			Y = Y_[i]
			Y_ptr = M * eBands[i]
		} else {
			Y = nil
		}
		N = M*eBands[i+1] - M*eBands[i]
		tell = int(ec.Tell_frac())

		if i != start {
			balance -= tell
		}
		remaining_bits = total_bits - tell - 1
		ctx.Remaining_bits = remaining_bits
		if i <= codedBands-1 {
			curr_balance = Inlines.Celt_sudiv(balance, Inlines.IMIN(3, codedBands-i))
			bits = Inlines.IMAX(0, Inlines.IMIN(16383, Inlines.IMIN(remaining_bits+1, pulses[i]+curr_balance)))
		} else {
			bits = 0
		}

		if resynth != 0 && M*eBands[i]-N >= M*eBands[start] && (update_lowband != 0 || lowband_offset == 0) {
			lowband_offset = i
		}

		ctx.Tf_change = tf_res[i]
		if i >= m.effEBands {
			X = norm
			X_ptr = 0
			if Y_ != nil {
				Y = norm
				Y_ptr = 0
			}
			lowband_scratch = nil
		}
		if i == end-1 {
			lowband_scratch = nil
		}

		if lowband_offset != 0 && (spread != Spread.SPREAD_AGGRESSIVE || B > 1 || tf_res[i] < 0) {
			fold_start := lowband_offset
			for M*eBands[fold_start] > effective_lowband+norm_offset {
				fold_start--
			}
			fold_end := lowband_offset - 1
			for M*eBands[fold_end] < effective_lowband+norm_offset+N {
				fold_end++
			}
			x_cm, y_cm := int64(0), int64(0)
			fold_i := fold_start
			for fold_i <= fold_end {
				x_cm |= int64(collapse_masks[fold_i*C+0])
				y_cm |= int64(collapse_masks[fold_i*C+C-1])
				fold_i++
			}
		} else {
			x_cm, y_cm := int64((1<<B)-1), int64((1<<B)-1)
		}

		if dual_stereo != 0 && i == intensity {
			if resynth != 0 {
				for j := 0; j < M*eBands[i]-norm_offset; j++ {
					norm[j] = Inlines.HALF32(norm[j] + norm[norm2+j])
				}
			}
			dual_stereo = 0
		}

		if dual_stereo != 0 {
			x_cm := int64(b.Quant_band(ctx, X, X_ptr, N, bits/2, B,
				effective_lowband != -1 ? norm : nil, effective_lowband, LM,
				last != 0 ? nil : norm, M*eBands[i]-norm_offset,
				CeltConstants.Q15ONE, lowband_scratch, lowband_scratch_ptr, int(x_cm)))
			y_cm := int64(b.Quant_band(ctx, Y, Y_ptr, N, bits/2, B,
				effective_lowband != -1 ? norm : nil, norm2+effective_lowband, LM,
				last != 0 ? nil : norm, norm2+(M*eBands[i]-norm_offset),
				CeltConstants.Q15ONE, lowband_scratch, lowband_scratch_ptr, int(y_cm)))
			collapse_masks[i*C+0] = int16(x_cm & 0xFF)
			collapse_masks[i*C+C-1] = int16(y_cm & 0xFF)
		} else {
			var cm_val int64
			if Y != nil {
				cm_val = int64(b.Quant_band_stereo(ctx, X, X_ptr, Y, Y_ptr, N, bits, B,
					effective_lowband != -1 ? norm : nil, effective_lowband, LM,
					last != 0 ? nil : norm, M*eBands[i]-norm_offset,
					lowband_scratch, lowband_scratch_ptr, int(x_cm|y_cm)))
			} else {
				cm_val = int64(b.Quant_band(ctx, X, X_ptr, N, bits, B,
					effective_lowband != -1 ? norm : nil, effective_lowband, LM,
					last != 0 ? nil : norm, M*eBands[i]-norm_offset,
					CeltConstants.Q15ONE, lowband_scratch, lowband_scratch_ptr, int(x_cm|y_cm)))
			}
			collapse_masks[i*C+0] = int16(cm_val & 0xFF)
			collapse_masks[i*C+C-1] = int16(cm_val & 0xFF)
		}
		balance += pulses[i] + tell
		update_lowband = Ternary(bits > (N<<EntropyCoder.BITRES), 1, 0)
	}
	seed.Val = ctx.Seed
}

// Helper function to simulate ternary operator
func Ternary(condition bool, trueVal, falseVal int) int {
	if condition {
		return trueVal
	}
	return falseVal
}