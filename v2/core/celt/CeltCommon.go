package celt

type CeltCommon struct {
	invTable []int16
}

func NewCeltCommon() *CeltCommon {
	return &CeltCommon{
		invTable: []int16{
			255, 255, 156, 110, 86, 70, 59, 51, 45, 40, 37, 33, 31, 28, 26, 25,
			23, 22, 21, 20, 19, 18, 17, 16, 16, 15, 15, 14, 13, 13, 12, 12,
			12, 12, 11, 11, 11, 10, 10, 10, 9, 9, 9, 9, 9, 9, 8, 8,
			8, 8, 8, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6,
			6, 6, 6, 6, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5,
			5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
			4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3,
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2,
		},
	}
}

func (c *CeltCommon) ComputeVbr(mode CeltMode, analysis AnalysisInfo, base_target int, LM int, bitrate int, lastCodedBands int, C int, intensity int,
	constrained_vbr int, stereo_saving int16, tot_boost int, tf_estimate int16, pitch_change int, maxDepth int16,
	variable_duration OpusFramesize, lfe int, has_surround_mask int, surround_masking int16, temporal_vbr int) int {

	target := base_target
	nbEBands := mode.nbEBands
	eBands := mode.eBands

	coded_bands := lastCodedBands
	if coded_bands == 0 {
		coded_bands = nbEBands
	}
	coded_bins := eBands[coded_bands] << LM
	if C == 2 {
		if intensity < coded_bands {
			coded_bins += eBands[intensity] << LM
		} else {
			coded_bins += eBands[coded_bands] << LM
		}
	}

	if analysis.enabled && analysis.valid != 0 && analysis.activity < .4 {
		boost := int(float32(coded_bins<<BITRES) * (.4 - analysis.activity))
		target -= boost
	}

	if C == 2 {
		coded_stereo_bands := c.min(intensity, coded_bands)
		coded_stereo_dof := (eBands[coded_stereo_bands] << LM) - coded_stereo_bands
		max_frac := (0.8 * float32(coded_stereo_dof)) / float32(coded_bins)
		if stereo_saving > 256 {
			stereo_saving = 256
		}

		stereo_boost := (int(stereo_saving)-26)*coded_stereo_dof<<BITRES >> 8
		if stereo_boost > int(max_frac*float32(target)) {
			stereo_boost = int(max_frac * float32(target))
		}
		target -= stereo_boost
	}

	target += tot_boost - (16 << LM)

	tf_calibration := int16(0)
	if variable_duration == OPUS_FRAMESIZE_VARIABLE {
		tf_calibration = 328 // QCONST16(0.02f, 14)
	} else {
		tf_calibration = 655 // QCONST16(0.04f, 14)
	}
	tf_boost := (int(tf_estimate)-int(tf_calibration))*target*2 >> 15
	target += tf_boost

	if analysis.enabled && analysis.valid != 0 && lfe == 0 {
		tonal := analysis.tonality
		if tonal < 0.15 {
			tonal = 0
		} else {
			tonal -= 0.15
		}
		tonal_boost := int(1.2 * float32(coded_bins<<BITRES) * (tonal - 0.09))
		if pitch_change != 0 {
			tonal_boost += int(0.8 * float32(coded_bins<<BITRES))
		}
		target += tonal_boost
	}

	if has_surround_mask != 0 && lfe == 0 {
		surround_boost := int(surround_masking) * (coded_bins << BITRES) >> DB_SHIFT
		if target/4 > surround_boost+base_target {
			target = target/4
		} else {
			target = base_target + surround_boost
		}
	}

	bins := eBands[nbEBands-2] << LM
	floor_depth := (C*bins<<BITRES)*int(maxDepth) >> DB_SHIFT
	if target>>2 > floor_depth {
		floor_depth = target >> 2
	}
	if floor_depth < target {
		target = floor_depth
	}

	if (has_surround_mask == 0 || lfe != 0) && (constrained_vbr != 0 || bitrate < 64000) {
		rate_factor := bitrate - 32000
		if rate_factor < 0 {
			rate_factor = 0
		}
		if constrained_vbr != 0 {
			if rate_factor > 21945 { // QCONST16(0.67f,15)
				rate_factor = 21945
			}
		}
		adj := rate_factor * (target - base_target) >> 15
		target = base_target + adj
	}

	if has_surround_mask == 0 && tf_estimate < 328 { // QCONST16(0.2f,14)
		amount := (c.max(0, c.min(32000, 96000-bitrate)) * 3) >> 15 // QCONST16(0.0000031f,30)
		tvbr_boost := amount * temporal_vbr * target >> 26
		target += tvbr_boost
	}

	if target > 2*base_target {
		target = 2 * base_target
	}

	return target
}

func (c *CeltCommon) TransientAnalysis(input [][]int, len int, C int, tf_estimate *int16, tf_chan *int) int {
	tmp := make([]int, len)
	mask_metric := 0
	*tf_chan = 0

	for ch := 0; ch < C; ch++ {
		var mem0, mem1 int
		unmask := 0
		for i := 0; i < len; i++ {
			x := input[ch][i] >> SIG_SHIFT
			y := x + mem0
			mem0 = mem1 + y - (x << 1)
			mem1 = x - (y >> 1)
			tmp[i] = y >> 2
		}
		for i := 0; i < 12; i++ {
			tmp[i] = 0
		}

		max := 0
		for i := 0; i < len; i++ {
			if c.abs(tmp[i]) > max {
				max = c.abs(tmp[i])
			}
		}
		shift := 14 - c.celtILog2(max+1)
		if shift != 0 {
			for i := 0; i < len; i++ {
				tmp[i] <<= shift
			}
		}

		sum := 0
		for i := 0; i < len; i++ {
			sum += (tmp[i] * tmp[i]) >> 16
		}

		mem0 = 0
		len2 := len / 2
		grouped := make([]int, len2)
		for i := 0; i < len2; i++ {
			val := ((tmp[2*i] * tmp[2*i]) + (tmp[2*i+1] * tmp[2*i+1])) >> 16
			grouped[i] = mem0 + ((val - mem0) >> 4)
			mem0 = grouped[i]
		}

		mem0 = 0
		maxE := 0
		for i := len2 - 1; i >= 0; i-- {
			grouped[i] = mem0 + ((grouped[i] - mem0) >> 3)
			mem0 = grouped[i]
			if grouped[i] > maxE {
				maxE = grouped[i]
			}
		}

		meanSqrt := c.celtSqrt(sum)
		maxESqrt := c.celtSqrt(maxE * len2 >> 1)
		norm := (len2 << 20) / (meanSqrt*maxESqrt>>1 + EPSILON)

		for i := 12; i < len2-5; i += 4 {
			id := grouped[i] * norm >> 15
			if id > 127 {
				id = 127
			}
			if id < 0 {
				id = 0
			}
			unmask += int(c.invTable[id])
		}
		unmask = unmask * 64 * 4 / (6 * (len2 - 17))
		if unmask > mask_metric {
			*tf_chan = ch
			mask_metric = unmask
		}
	}

	is_transient := 0
	if mask_metric > 200 {
		is_transient = 1
	}

	tf_max := 0
	if mask_metric > 0 {
		tf_max = (c.celtSqrt(27*mask_metric) - 42)
	}
	if tf_max < 0 {
		tf_max = 0
	}

	constVal := int16(113) // QCONST16(0.0069f,14)
	*tf_estimate = int16(c.celtSqrt(int(c.max(0, (int(constVal)*tf_max<<14)-233203)))
	return is_transient
}

func (c *CeltCommon) PatchTransientDecision(newE [][]int, oldE [][]int, nbEBands int, start int, end int, C int) int {
	mean_diff := 0
	spread_old := make([]int, 26)

	if C == 1 {
		spread_old[start] = oldE[0][start]
		for i := start + 1; i < end; i++ {
			val := spread_old[i-1] - (1 << DB_SHIFT)
			spread_old[i] = c.max(val, oldE[0][i])
		}
	} else {
		spread_old[start] = c.max(oldE[0][start], oldE[1][start])
		for i := start + 1; i < end; i++ {
			val := spread_old[i-1] - (1 << DB_SHIFT)
			maxOld := c.max(oldE[0][i], oldE[1][i])
			spread_old[i] = c.max(val, maxOld)
		}
	}

	for i := end - 2; i >= start; i-- {
		val := spread_old[i+1] - (1 << DB_SHIFT)
		spread_old[i] = c.max(spread_old[i], val)
	}

	startIdx := c.max(2, start)

	for ch := 0; ch < C; ch++ {
		for i := startIdx; i < end-1; i++ {
			x1 := c.max(0, newE[ch][i])
			x2 := c.max(0, spread_old[i])
			diff := x1 - x2
			if diff > 0 {
				mean_diff += diff
			}
		}
	}
	mean_diff /= C * (end - 1 - startIdx)

	if mean_diff > (1 << DB_SHIFT) {
		return 1
	}
	return 0
}

// Helper methods
func (c *CeltCommon) abs(x int) int {
	if x < 0 {
		return -x
	}
	return x
}

func (c *CeltCommon) celtSqrt(x int) int {
	return int(c.iSqrt32(uint32(x<<10))) >> 5
}

func (c *CeltCommon) celtILog2(x int) int {
	if x <= 0 {
		return 0
	}
	log := 0
	for x > 1 {
		log++
		x >>= 1
	}
	return log
}

func (c *CeltCommon) min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func (c *CeltCommon) max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func (c *CeltCommon) iSqrt32(val uint32) uint32 {
	if val == 0 {
		return 0
	}
	g := uint32(0x8000)
	for i := 0; i < 16; i++ {
		temp := (0x8000 | g) >> 1
		if temp*temp <= val>>(2*i) {
			g = temp
		} else {
			g >>= 1
		}
	}
	return g
}
 /**
     * Apply window and compute the MDCT for all sub-frames and all channels in
     * a frame
     */
    static void compute_mdcts(CeltMode mode, int shortBlocks, int[][] input,
            int[][] output, int C, int CC, int LM, int upsample) {
        int overlap = mode.overlap;
        int N;
        int B;
        int shift;
        int i, b, c;
        if (shortBlocks != 0) {
            B = shortBlocks;
            N = mode.shortMdctSize;
            shift = mode.maxLM;
        } else {
            B = 1;
            N = mode.shortMdctSize << LM;
            shift = mode.maxLM - LM;
        }
        c = 0;
        do {
            for (b = 0; b < B; b++) {
                /* Interleaving the sub-frames while doing the MDCTs */
                MDCT.clt_mdct_forward(
                        mode.mdct,
                        input[c],
                        b * N,
                        output[c],
                        b,
                        mode.window,
                        overlap,
                        shift,
                        B);
            }
        } while (++c < CC);

        if (CC == 2 && C == 1) {
            for (i = 0; i < B * N; i++) {
                output[0][i] = Inlines.ADD32(Inlines.HALF32(output[0][i]), Inlines.HALF32(output[1][i]));
            }
        }
        if (upsample != 1) {
            c = 0;
            do {
                int bound = B * N / upsample;
                for (i = 0; i < bound; i++) {
                    output[c][i] *= upsample;
                }
                Arrays.MemSetWithOffset(output[c], 0, bound, B * N - bound);
            } while (++c < C);
        }
    }

    static void celt_preemphasis(short[] pcmp, int pcmp_ptr, int[] inp, int inp_ptr,
            int N, int CC, int upsample, int[] coef, BoxedValueInt mem, int clip) {
        int i;
        int coef0;
        int m;
        int Nu;

        coef0 = coef[0];
        m = mem.Val;

        /* Fast path for the normal 48kHz case and no clipping */
        if (coef[1] == 0 && upsample == 1 && clip == 0) {
            for (i = 0; i < N; i++) {
                int x = pcmp[pcmp_ptr + (CC * i)];
                /* Apply pre-emphasis */
                inp[inp_ptr + i] = Inlines.SHL32(x, CeltConstants.SIG_SHIFT) - m;
                m = Inlines.SHR32(Inlines.MULT16_16(coef0, x), 15 - CeltConstants.SIG_SHIFT);
            }
            mem.Val = m;
            return;
        }

        Nu = N / upsample;
        if (upsample != 1) {
            Arrays.MemSetWithOffset(inp, 0, inp_ptr, N);
        }
        for (i = 0; i < Nu; i++) {
            inp[inp_ptr + (i * upsample)] = pcmp[pcmp_ptr + (CC * i)];
        }

        for (i = 0; i < N; i++) {
            int x;
            x = (inp[inp_ptr + i]);
            /* Apply pre-emphasis */
            inp[inp_ptr + i] = Inlines.SHL32(x, CeltConstants.SIG_SHIFT) - m;
            m = Inlines.SHR32(Inlines.MULT16_16(coef0, x), 15 - CeltConstants.SIG_SHIFT);
        }

        mem.Val = m;
    }

    static void celt_preemphasis(short[] pcmp, int[] inp, int inp_ptr,
            int N, int CC, int upsample, int[] coef, BoxedValueInt mem, int clip) {
        int i;
        int coef0;
        int m;
        int Nu;

        coef0 = coef[0];
        m = mem.Val;

        /* Fast path for the normal 48kHz case and no clipping */
        if (coef[1] == 0 && upsample == 1 && clip == 0) {
            for (i = 0; i < N; i++) {
                int x;
                x = pcmp[CC * i];
                /* Apply pre-emphasis */
                inp[inp_ptr + i] = Inlines.SHL32(x, CeltConstants.SIG_SHIFT) - m;
                m = Inlines.SHR32(Inlines.MULT16_16(coef0, x), 15 - CeltConstants.SIG_SHIFT);
            }
            mem.Val = m;
            return;
        }

        Nu = N / upsample;
        if (upsample != 1) {
            Arrays.MemSetWithOffset(inp, 0, inp_ptr, N);
        }
        for (i = 0; i < Nu; i++) {
            inp[inp_ptr + (i * upsample)] = pcmp[CC * i];
        }

        for (i = 0; i < N; i++) {
            int x;
            x = (inp[inp_ptr + i]);
            /* Apply pre-emphasis */
            inp[inp_ptr + i] = Inlines.SHL32(x, CeltConstants.SIG_SHIFT) - m;
            m = Inlines.SHR32(Inlines.MULT16_16(coef0, x), 15 - CeltConstants.SIG_SHIFT);
        }

        mem.Val = m;
    }

    static int l1_metric(int[] tmp, int N, int LM, int bias) {
        int i;
        int L1;
        L1 = 0;
        for (i = 0; i < N; i++) {
            L1 += Inlines.EXTEND32(Inlines.ABS32(tmp[i]));
        }

        /* When in doubt, prefer good freq resolution */
        L1 = Inlines.MAC16_32_Q15(L1, (LM * bias), (L1));
        return L1;

    }

    static int tf_analysis(CeltMode m, int len, int isTransient,
            int[] tf_res, int lambda, int[][] X, int N0, int LM,
            BoxedValueInt tf_sum, int tf_estimate, int tf_chan) {
        int i;
        int[] metric;
        int cost0;
        int cost1;
        int[] path0;
        int[] path1;
        int[] tmp;
        int[] tmp_1;
        int sel;
        int[] selcost = new int[2];
        int tf_select = 0;
        int bias;

        bias = Inlines.MULT16_16_Q14(((short) (0.5 + (.04f) * (((int) 1) << (15))))/*Inlines.QCONST16(.04f, 15)*/, Inlines.MAX16((short) (0 - ((short) (0.5 + (.25f) * (((int) 1) << (14))))/*Inlines.QCONST16(.25f, 14)*/), (((short) (0.5 + (.5f) * (((int) 1) << (14))))/*Inlines.QCONST16(.5f, 14)*/ - tf_estimate)));
        /*printf("%f ", bias);*/

        metric = new int[len];
        tmp = new int[(m.eBands[len] - m.eBands[len - 1]) << LM];
        tmp_1 = new int[(m.eBands[len] - m.eBands[len - 1]) << LM];
        path0 = new int[len];
        path1 = new int[len];

        tf_sum.Val = 0;
        for (i = 0; i < len; i++) {
            int k, N;
            int narrow;
            int L1, best_L1;
            int best_level = 0;
            N = (m.eBands[i + 1] - m.eBands[i]) << LM;
            /* band is too narrow to be split down to LM=-1 */
            narrow = ((m.eBands[i + 1] - m.eBands[i]) == 1) ? 1 : 0;
            System.arraycopy(X[tf_chan], (m.eBands[i] << LM), tmp, 0, N);
            /* Just add the right channel if we're in stereo */
 /*if (C==2)
               for (j=0;j<N;j++)
                  tmp[j] = ADD16(SHR16(tmp[j], 1),SHR16(X[N0+j+(m.eBands[i]<<LM)], 1));*/
            L1 = l1_metric(tmp, N, isTransient != 0 ? LM : 0, bias);
            best_L1 = L1;
            /* Check the -1 case for transients */
            if (isTransient != 0 && narrow == 0) {
                System.arraycopy(tmp, 0, tmp_1, 0, N);
                Bands.haar1ZeroOffset(tmp_1, N >> LM, 1 << LM);
                L1 = l1_metric(tmp_1, N, LM + 1, bias);
                if (L1 < best_L1) {
                    best_L1 = L1;
                    best_level = -1;
                }
            }
            /*printf ("%f ", L1);*/
            for (k = 0; k < LM + (!(isTransient != 0 || narrow != 0) ? 1 : 0); k++) {
                int B;

                if (isTransient != 0) {
                    B = (LM - k - 1);
                } else {
                    B = k + 1;
                }

                Bands.haar1ZeroOffset(tmp, N >> k, 1 << k);

                L1 = l1_metric(tmp, N, B, bias);

                if (L1 < best_L1) {
                    best_L1 = L1;
                    best_level = k + 1;
                }
            }
            /*printf ("%d ", isTransient ? LM-best_level : best_level);*/
 /* metric is in Q1 to be able to select the mid-point (-0.5) for narrower bands */
            if (isTransient != 0) {
                metric[i] = 2 * best_level;
            } else {
                metric[i] = -2 * best_level;
            }
            tf_sum.Val += (isTransient != 0 ? LM : 0) - metric[i] / 2;
            /* For bands that can't be split to -1, set the metric to the half-way point to avoid
               biasing the decision */
            if (narrow != 0 && (metric[i] == 0 || metric[i] == -2 * LM)) {
                metric[i] -= 1;
            }
            /*printf("%d ", metric[i]);*/
        }
        /*printf("\n");*/
 /* Search for the optimal tf resolution, including tf_select */
        tf_select = 0;
        for (sel = 0; sel < 2; sel++) {
            cost0 = 0;
            cost1 = isTransient != 0 ? 0 : lambda;
            for (i = 1; i < len; i++) {
                int curr0, curr1;
                curr0 = Inlines.IMIN(cost0, cost1 + lambda);
                curr1 = Inlines.IMIN(cost0 + lambda, cost1);
                cost0 = curr0 + Inlines.abs(metric[i] - 2 * CeltTables.tf_select_table[LM][4 * isTransient + 2 * sel + 0]);
                cost1 = curr1 + Inlines.abs(metric[i] - 2 * CeltTables.tf_select_table[LM][4 * isTransient + 2 * sel + 1]);
            }
            cost0 = Inlines.IMIN(cost0, cost1);
            selcost[sel] = cost0;
        }
        /* For now, we're conservative and only allow tf_select=1 for transients.
         * If tests confirm it's useful for non-transients, we could allow it. */
        if (selcost[1] < selcost[0] && isTransient != 0) {
            tf_select = 1;
        }
        cost0 = 0;
        cost1 = isTransient != 0 ? 0 : lambda;
        /* Viterbi forward pass */
        for (i = 1; i < len; i++) {
            int curr0, curr1;
            int from0, from1;

            from0 = cost0;
            from1 = cost1 + lambda;
            if (from0 < from1) {
                curr0 = from0;
                path0[i] = 0;
            } else {
                curr0 = from1;
                path0[i] = 1;
            }

            from0 = cost0 + lambda;
            from1 = cost1;
            if (from0 < from1) {
                curr1 = from0;
                path1[i] = 0;
            } else {
                curr1 = from1;
                path1[i] = 1;
            }
            cost0 = curr0 + Inlines.abs(metric[i] - 2 * CeltTables.tf_select_table[LM][4 * isTransient + 2 * tf_select + 0]);
            cost1 = curr1 + Inlines.abs(metric[i] - 2 * CeltTables.tf_select_table[LM][4 * isTransient + 2 * tf_select + 1]);
        }
        tf_res[len - 1] = cost0 < cost1 ? 0 : 1;
        /* Viterbi backward pass to check the decisions */
        for (i = len - 2; i >= 0; i--) {
            if (tf_res[i + 1] == 1) {
                tf_res[i] = path1[i + 1];
            } else {
                tf_res[i] = path0[i + 1];
            }
        }
        /*printf("%d %f\n", *tf_sum, tf_estimate);*/

        return tf_select;
    }

    static void tf_encode(int start, int end, int isTransient, int[] tf_res, int LM, int tf_select, EntropyCoder enc) {
        int curr, i;
        int tf_select_rsv;
        int tf_changed;
        int logp;
        int budget;
        int tell;
        budget = enc.storage * 8;
        tell = enc.tell();
        logp = isTransient != 0 ? 2 : 4;
        /* Reserve space to code the tf_select decision. */
        tf_select_rsv = (LM > 0 && tell + logp + 1 <= budget) ? 1 : 0;
        budget -= tf_select_rsv;
        curr = tf_changed = 0;
        for (i = start; i < end; i++) {
            if (tell + logp <= budget) {
                enc.enc_bit_logp(tf_res[i] ^ curr, logp);
                tell = enc.tell();
                curr = tf_res[i];
                tf_changed |= curr;
            } else {
                tf_res[i] = curr;
            }
            logp = isTransient != 0 ? 4 : 5;
        }
        /* Only code tf_select if it would actually make a difference. */
        if (tf_select_rsv != 0
                && CeltTables.tf_select_table[LM][4 * isTransient + 0 + tf_changed]
                != CeltTables.tf_select_table[LM][4 * isTransient + 2 + tf_changed]) {
            enc.enc_bit_logp(tf_select, 1);
        } else {
            tf_select = 0;
        }
        for (i = start; i < end; i++) {
            tf_res[i] = CeltTables.tf_select_table[LM][4 * isTransient + 2 * tf_select + tf_res[i]];
        }
        /*for(i=0;i<end;i++)printf("%d ", isTransient ? tf_res[i] : LM+tf_res[i]);printf("\n");*/
    }

    static int alloc_trim_analysis(CeltMode m, int[][] X,
            int[][] bandLogE, int end, int LM, int C,
            AnalysisInfo analysis, BoxedValueInt stereo_saving, int tf_estimate,
            int intensity, int surround_trim) {
        int i;
        int diff = 0;
        int c;
        int trim_index;
        int trim = ((short) (0.5 + (5.0f) * (((int) 1) << (8))))/*Inlines.QCONST16(5.0f, 8)*/;
        int logXC, logXC2;
        if (C == 2) {
            int sum = 0;
            /* Q10 */
            int minXC;
            /* Q10 */
 /* Compute inter-channel correlation for low frequencies */
            for (i = 0; i < 8; i++) {
                int partial;
                partial = Kernels.celt_inner_prod(X[0], (m.eBands[i] << LM), X[1], (m.eBands[i] << LM),
                        (m.eBands[i + 1] - m.eBands[i]) << LM);
                sum = Inlines.ADD16(sum, Inlines.EXTRACT16(Inlines.SHR32(partial, 18)));
            }
            sum = Inlines.MULT16_16_Q15(((short) (0.5 + (1.0f / 8) * (((int) 1) << (15))))/*Inlines.QCONST16(1.0f / 8, 15)*/, sum);
            sum = Inlines.MIN16(((short) (0.5 + (1.0f) * (((int) 1) << (10))))/*Inlines.QCONST16(1.0f, 10)*/, Inlines.ABS32(sum));
            minXC = sum;
            for (i = 8; i < intensity; i++) {
                int partial;
                partial = Kernels.celt_inner_prod(X[0], (m.eBands[i] << LM), X[1], (m.eBands[i] << LM),
                        (m.eBands[i + 1] - m.eBands[i]) << LM);
                minXC = Inlines.MIN16(minXC, Inlines.ABS16(Inlines.EXTRACT16(Inlines.SHR32(partial, 18))));
            }
            minXC = Inlines.MIN16(((short) (0.5 + (1.0f) * (((int) 1) << (10))))/*Inlines.QCONST16(1.0f, 10)*/, Inlines.ABS32(minXC));
            /*printf ("%f\n", sum);*/
 /* mid-side savings estimations based on the LF average*/
            logXC = Inlines.celt_log2(((int) (0.5 + (1.001f) * (((int) 1) << (20))))/*Inlines.QCONST32(1.001f, 20)*/ - Inlines.MULT16_16(sum, sum));
            /* mid-side savings estimations based on min correlation */
            logXC2 = Inlines.MAX16(Inlines.HALF16(logXC), Inlines.celt_log2(((int) (0.5 + (1.001f) * (((int) 1) << (20))))/*Inlines.QCONST32(1.001f, 20)*/ - Inlines.MULT16_16(minXC, minXC)));
            /* Compensate for Q20 vs Q14 input and convert output to Q8 */
            logXC = (Inlines.PSHR32(logXC - ((short) (0.5 + (6.0f) * (((int) 1) << (CeltConstants.DB_SHIFT))))/*Inlines.QCONST16(6.0f, CeltConstants.DB_SHIFT)*/, CeltConstants.DB_SHIFT - 8));
            logXC2 = (Inlines.PSHR32(logXC2 - ((short) (0.5 + (6.0f) * (((int) 1) << (CeltConstants.DB_SHIFT))))/*Inlines.QCONST16(6.0f, CeltConstants.DB_SHIFT)*/, CeltConstants.DB_SHIFT - 8));

            trim += Inlines.MAX16((0 - ((short) (0.5 + (4.0f) * (((int) 1) << (8))))/*Inlines.QCONST16(4.0f, 8)*/), Inlines.MULT16_16_Q15(((short) (0.5 + (.75f) * (((int) 1) << (15))))/*Inlines.QCONST16(.75f, 15)*/, logXC));
            stereo_saving.Val = Inlines.MIN16((stereo_saving.Val + ((short) (0.5 + (0.25f) * (((int) 1) << (8))))/*Inlines.QCONST16(0.25f, 8)*/), (0 - Inlines.HALF16(logXC2)));
        }

        /* Estimate spectral tilt */
        c = 0;
        do {
            for (i = 0; i < end - 1; i++) {
                diff += bandLogE[c][i] * (int) (2 + 2 * i - end);
            }
        } while (++c < C);
        diff /= C * (end - 1);
        /*printf("%f\n", diff);*/
        trim -= Inlines.MAX16(Inlines.NEG16(((short) (0.5 + (2.0f) * (((int) 1) << (8))))/*Inlines.QCONST16(2.0f, 8)*/), Inlines.MIN16(((short) (0.5 + (2.0f) * (((int) 1) << (8))))/*Inlines.QCONST16(2.0f, 8)*/, (Inlines.SHR16((diff + ((short) (0.5 + (1.0f) * (((int) 1) << (CeltConstants.DB_SHIFT))))/*Inlines.QCONST16(1.0f, CeltConstants.DB_SHIFT)*/), CeltConstants.DB_SHIFT - 8) / 6)));
        trim -= Inlines.SHR16(surround_trim, CeltConstants.DB_SHIFT - 8);
        trim = (trim - 2 * Inlines.SHR16(tf_estimate, 14 - 8));
        if (analysis.enabled && analysis.valid != 0) {
            trim -= Inlines.MAX16(-((short) (0.5 + (2.0f) * (((int) 1) << (8))))/*Inlines.QCONST16(2.0f, 8)*/, Inlines.MIN16(((short) (0.5 + (2.0f) * (((int) 1) << (8))))/*Inlines.QCONST16(2.0f, 8)*/,
                            (int) (((short) (0.5 + (2.0f) * (((int) 1) << (8))))/*Inlines.QCONST16(2.0f, 8)*/ * (analysis.tonality_slope + .05f))));
        }
        trim_index = Inlines.PSHR32(trim, 8);
        trim_index = Inlines.IMAX(0, Inlines.IMIN(10, trim_index));
        /*printf("%d\n", trim_index);*/

        return trim_index;
    }

    static int stereo_analysis(CeltMode m, int[][] X,
            int LM) {
        int i;
        int thetas;
        int sumLR = CeltConstants.EPSILON, sumMS = CeltConstants.EPSILON;

        /* Use the L1 norm to model the entropy of the L/R signal vs the M/S signal */
        for (i = 0; i < 13; i++) {
            int j;
            for (j = m.eBands[i] << LM; j < m.eBands[i + 1] << LM; j++) {
                int L, R, M, S;
                /* We cast to 32-bit first because of the -32768 case */
                L = Inlines.EXTEND32(X[0][j]);
                R = Inlines.EXTEND32(X[1][j]);
                M = Inlines.ADD32(L, R);
                S = Inlines.SUB32(L, R);
                sumLR = Inlines.ADD32(sumLR, Inlines.ADD32(Inlines.ABS32(L), Inlines.ABS32(R)));
                sumMS = Inlines.ADD32(sumMS, Inlines.ADD32(Inlines.ABS32(M), Inlines.ABS32(S)));
            }
        }
        sumMS = Inlines.MULT16_32_Q15(((short) (0.5 + (0.707107f) * (((int) 1) << (15))))/*Inlines.QCONST16(0.707107f, 15)*/, sumMS);
        thetas = 13;
        /* We don't need thetas for lower bands with LM<=1 */
        if (LM <= 1) {
            thetas -= 8;
        }
        return (Inlines.MULT16_32_Q15(((m.eBands[13] << (LM + 1)) + thetas), sumMS)
                > Inlines.MULT16_32_Q15((m.eBands[13] << (LM + 1)), sumLR)) ? 1 : 0;
    }
func (c *CeltCommon) MedianOf5(x []int, offset int) int {
	t0, t1, t2, t3, t4 := 0, 0, 0, 0, 0
	t2 = x[offset+2]
	
	if x[offset] > x[offset+1] {
		t0 = x[offset+1]
		t1 = x[offset]
	} else {
		t0 = x[offset]
		t1 = x[offset+1]
	}
	
	if x[offset+3] > x[offset+4] {
		t3 = x[offset+4]
		t4 = x[offset+3]
	} else {
		t3 = x[offset+3]
		t4 = x[offset+4]
	}
	
	if t0 > t3 {
		t3, t0 = t0, t3
		t4, t1 = t1, t4
	}
	
	if t2 > t1 {
		if t1 < t3 {
			return c.min(t2, t3)
		}
		return c.min(t4, t1)
	} else if t2 < t3 {
		return c.min(t1, t3)
	}
	return c.min(t2, t4)
}

func (c *CeltCommon) MedianOf3(x []int, offset int) int {
	t0, t1 := 0, 0
	if x[offset] > x[offset+1] {
		t0 = x[offset+1]
		t1 = x[offset]
	} else {
		t0 = x[offset]
		t1 = x[offset+1]
	}
	t2 := x[offset+2]
	
	if t1 < t2 {
		return t1
	} else if t0 < t2 {
		return t2
	}
	return t0
}

func (c *CeltCommon) DynallocAnalysis(bandLogE [][]int, bandLogE2 [][]int, nbEBands int, start int, end int, channels int, offsets []int, lsbDepth int, logN []int16, isTransient int, vbr int, constrainedVbr int, eBands []int, LM int, effectiveBytes int, totBoost *int, lfe int, surroundDynalloc []int) int {
	tot := 0
	maxDepth := int(-31.9 * float32(1<<CeltConstants.DB_SHIFT)) // QCONST16(31.9f, DB_SHIFT)
	noiseFloor := make([]int, channels*nbEBands)
	follower := make([][]int, channels)
	for i := range follower {
		follower[i] = make([]int, nbEBands)
	}

	for i := 0; i < end; i++ {
		noiseFloor[i] = int(logN[i])*int(0.0625*float32(1<<CeltConstants.DB_SHIFT))>>16 +
			int(0.5*float32(1<<CeltConstants.DB_SHIFT)) +
			(9-lsbDepth)<<CeltConstants.DB_SHIFT -
			int(CeltTables.EMeans[i])<<6 +
			int(0.0062*float32(1<<CeltConstants.DB_SHIFT))*((i+5)*(i+5))>>16
	}

	for ch := 0; ch < channels; ch++ {
		for i := 0; i < end; i++ {
			if bandLogE[ch][i]-noiseFloor[i] > maxDepth {
				maxDepth = bandLogE[ch][i] - noiseFloor[i]
			}
		}
	}

	if effectiveBytes > 50 && LM >= 1 && lfe == 0 {
		last := 0
		for ch := 0; ch < channels; ch++ {
			f := follower[ch]
			f[0] = bandLogE2[ch][0]
			for i := 1; i < end; i++ {
				if bandLogE2[ch][i] > bandLogE2[ch][i-1]+int(0.5*float32(1<<CeltConstants.DB_SHIFT)) {
					last = i
				}
				f[i] = c.min(f[i-1]+int(1.5*float32(1<<CeltConstants.DB_SHIFT)), bandLogE2[ch][i])
			}
			for i := last - 1; i >= 0; i-- {
				f[i] = c.min(f[i], c.min(f[i+1]+int(2.0*float32(1<<CeltConstants.DB_SHIFT)), bandLogE2[ch][i]))
			}

			offsetVal := int(1.0 * float32(1<<CeltConstants.DB_SHIFT))
			for i := 2; i < end-2; i++ {
				med := c.MedianOf5(bandLogE2[ch], i-2)
				f[i] = c.max(f[i], med-offsetVal)
			}
			med := c.MedianOf3(bandLogE2[ch], 0)
			f[0] = c.max(f[0], med-offsetVal)
			f[1] = c.max(f[1], med-offsetVal)
			med = c.MedianOf3(bandLogE2[ch], end-3)
			f[end-2] = c.max(f[end-2], med-offsetVal)
			f[end-1] = c.max(f[end-1], med-offsetVal)

			for i := 0; i < end; i++ {
				f[i] = c.max(f[i], noiseFloor[i])
			}
		}

		if channels == 2 {
			for i := start; i < end; i++ {
				follower[1][i] = c.max(follower[1][i], follower[0][i]-int(4.0*float32(1<<CeltConstants.DB_SHIFT)))
				follower[0][i] = c.max(follower[0][i], follower[1][i]-int(4.0*float32(1<<CeltConstants.DB_SHIFT)))
				follower[0][i] = (c.max(0, bandLogE[0][i]-follower[0][i]) + c.max(0, bandLogE[1][i]-follower[1][i])) / 2
			}
		} else {
			for i := start; i < end; i++ {
				follower[0][i] = c.max(0, bandLogE[0][i]-follower[0][i])
			}
		}

		for i := start; i < end; i++ {
			follower[0][i] = c.max(follower[0][i], surroundDynalloc[i])
		}

		if (vbr == 0 || constrainedVbr != 0) && isTransient == 0 {
			for i := start; i < end; i++ {
				follower[0][i] /= 2
			}
		}

		for i := start; i < end; i++ {
			width := channels * (eBands[i+1] - eBands[i]) << LM
			boost := follower[0][i] >> CeltConstants.DB_SHIFT
			var boostBits int

			switch {
			case width < 6:
				boostBits = boost * width << EntropyCoder.BITRES
			case width > 48:
				boostBits = (boost * width << EntropyCoder.BITRES) / 8
			default:
				boostBits = boost * 6 << EntropyCoder.BITRES
			}

			if (vbr == 0 || (constrainedVbr != 0 && isTransient == 0)) && 
				(tot+boostBits)>>EntropyCoder.BITRES>>3 > effectiveBytes/4 {
				cap := (effectiveBytes / 4) << EntropyCoder.BITRES << 3
				offsets[i] = cap - tot
				tot = cap
				break
			} else {
				offsets[i] = boost
				tot += boostBits
			}
		}
	}

	*totBoost = tot
	return maxDepth
}

func (c *CeltCommon) Deemphasis(input [][]int, inputPtrs []int, pcm []int16, pcmOffset int, N int, channels int, downsample int, coef []int16, mem []int, accum int) {
	Nd := N / downsample
	coef0 := int(coef[0])
	scratch := make([]int, N)

	for ch := 0; ch < channels; ch++ {
		m := mem[ch]
		x := input[ch]
		xPtr := inputPtrs[ch]
		yOffset := pcmOffset + ch

		if downsample > 1 {
			for j := 0; j < N; j++ {
				tmp := x[xPtr+j] + m + CeltConstants.VERY_SMALL
				m = int((int32(coef0) * int32(tmp)) >> 15)
				scratch[j] = tmp
			}
		} else if accum != 0 {
			for j := 0; j < N; j++ {
				tmp := x[xPtr+j] + m + CeltConstants.VERY_SMALL
				m = int((int32(coef0) * int32(tmp)) >> 15)
				idx := yOffset + j*channels
				pcm[idx] = int16(c.SAT16(int(pcm[idx]) + (tmp >> CeltConstants.SIG_SHIFT)))
			}
		} else {
			for j := 0; j < N; j++ {
				tmp := x[xPtr+j] + m + CeltConstants.VERY_SMALL
				if x[xPtr+j] > 0 && m > 0 && tmp < 0 {
					tmp = math.MaxInt32
					m = math.MaxInt32
				} else {
					m = int((int32(coef0) * int32(tmp)) >> 15)
				}
				pcm[yOffset+j*channels] = int16(tmp >> CeltConstants.SIG_SHIFT)
			}
		}
		mem[ch] = m

		if downsample > 1 {
			for j := 0; j < Nd; j++ {
				pcm[yOffset+j*channels] = int16(scratch[j*downsample] >> CeltConstants.SIG_SHIFT)
			}
		}
	}
}

func (c *CeltCommon) CeltSynthesis(mode *CeltMode, X [][]int, outSyn [][]int, outSynPtrs []int, oldBandE []int, start int, effEnd int, channels int, CC int, isTransient int, LM int, downsample int, silence int) {
	overlap := mode.overlap
	nbEBands := mode.nbEBands
	N := mode.shortMdctSize << LM
	freq := make([]int, N)
	M := 1 << LM

	var B, NB, shift int
	if isTransient != 0 {
		B = M
		NB = mode.shortMdctSize
		shift = mode.maxLM
	} else {
		B = 1
		NB = mode.shortMdctSize << LM
		shift = mode.maxLM - LM
	}

	if CC == 2 && channels == 1 {
		c.Bands.DenormaliseBands(mode, X[0], freq, 0, oldBandE, 0, start, effEnd, M, downsample, silence)
		freq2 := outSynPtrs[1] + overlap/2
		copy(outSyn[1][freq2:freq2+N], freq[:N])
		
		for b := 0; b < B; b++ {
			c.mdct.CltMdctBackward(outSyn[1], freq2+b, outSyn[0], outSynPtrs[0]+NB*b, mode.window, overlap, shift, B)
		}
		for b := 0; b < B; b++ {
			c.mdct.CltMdctBackward(freq, b, outSyn[1], outSynPtrs[1]+NB*b, mode.window, overlap, shift, B)
		}
	} else if CC == 1 && channels == 2 {
		freq2 := outSynPtrs[0] + overlap/2
		c.Bands.DenormaliseBands(mode, X[0], freq, 0, oldBandE, 0, start, effEnd, M, downsample, silence)
		c.Bands.DenormaliseBands(mode, X[1], outSyn[0], freq2, oldBandE, nbEBands, start, effEnd, M, downsample, silence)
		
		for i := 0; i < N; i++ {
			freq[i] = (freq[i] + outSyn[0][freq2+i]) / 2
		}
		for b := 0; b < B; b++ {
			c.mdct.CltMdctBackward(freq, b, outSyn[0], outSynPtrs[0]+NB*b, mode.window, overlap, shift, B)
		}
	} else {
		for ch := 0; ch < CC; ch++ {
			c.Bands.DenormaliseBands(mode, X[ch], freq, 0, oldBandE, ch*nbEBands, start, effEnd, M, downsample, silence)
			for b := 0; b < B; b++ {
				c.mdct.CltMdctBackward(freq, b, outSyn[ch], outSynPtrs[ch]+NB*b, mode.window, overlap, shift, B)
			}
		}
	}
}

// Helper function for saturation
func (c *CeltCommon) SAT16(x int) int {
	if x < -32768 {
		return -32768
	}
	if x > 32767 {
		return 32767
	}
	return x
}
func (c *CeltCommon) TfDecode(start int, end int, isTransient int, tfRes []int, LM int, dec *EntropyCoder) {
	budget := dec.Storage * 8
	tell := dec.Tell()
	logp := 2
	if isTransient == 0 {
		logp = 4
	}
	
	tfSelectRsv := 0
	if LM > 0 && tell+logp+1 <= budget {
		tfSelectRsv = 1
	}
	budget -= tfSelectRsv
	
	curr := 0
	tfChanged := 0
	for i := start; i < end; i++ {
		if tell+logp <= budget {
			bit := dec.DecBitLogp(logp)
			curr ^= bit
			tell = dec.Tell()
			if bit != 0 {
				tfChanged = 1
			}
		}
		tfRes[i] = curr
		if isTransient != 0 {
			logp = 4
		} else {
			logp = 5
		}
	}
	
	tfSelect := 0
	if tfSelectRsv != 0 {
		idx1 := 4*isTransient + 0 + tfChanged
		idx2 := 4*isTransient + 2 + tfChanged
		if c.tfSelectTable[LM][idx1] != c.tfSelectTable[LM][idx2] {
			tfSelect = dec.DecBitLogp(1)
		}
	}
	
	for i := start; i < end; i++ {
		tfRes[i] = c.tfSelectTable[LM][4*isTransient+2*tfSelect+tfRes[i]]
	}
}

func (c *CeltCommon) CeltPlcPitchSearch(decodeMem [][]int, channels int) int {
	lpPitchBuf := make([]int, CeltConstants.DECODE_BUFFER_SIZE>>1)
	c.Pitch.PitchDownsample(decodeMem, lpPitchBuf, CeltConstants.DECODE_BUFFER_SIZE, channels)
	
	pitchIndex := 0
	c.Pitch.PitchSearch(
		lpPitchBuf,
		CeltConstants.PLC_PITCH_LAG_MAX>>1,
		lpPitchBuf,
		CeltConstants.DECODE_BUFFER_SIZE-CeltConstants.PLC_PITCH_LAG_MAX,
		CeltConstants.PLC_PITCH_LAG_MAX-CeltConstants.PLC_PITCH_LAG_MIN,
		&pitchIndex,
	)
	
	return CeltConstants.PLC_PITCH_LAG_MAX - pitchIndex
}

func (c *CeltCommon) ResamplingFactor(rate int) int {
	switch rate {
	case 48000:
		return 1
	case 24000:
		return 2
	case 16000:
		return 3
	case 12000:
		return 4
	case 8000:
		return 6
	default:
		panic("unsupported sample rate")
	}
}

func (c *CeltCommon) CombFilterConst(y []int, yPtr int, x []int, xPtr int, T int, N int, g10 int, g11 int, g12 int) {
	xpt := xPtr - T
	x4 := x[xpt-2]
	x3 := x[xpt-1]
	x2 := x[xpt]
	x1 := x[xpt+1]
	
	for i := 0; i < N; i++ {
		x0 := x[xpt+i+2]
		y[yPtr+i] = x[xPtr+i] +
			c.Mult16_32_Q15(g10, x2) +
			c.Mult16_32_Q15(g11, x1+x3) +
			c.Mult16_32_Q15(g12, x0+x4)
		
		x4 = x3
		x3 = x2
		x2 = x1
		x1 = x0
	}
}

func (c *CeltCommon) CombFilter(y []int, yPtr int, x []int, xPtr int, T0 int, T1 int, N int,
	g0 int, g1 int, tapset0 int, tapset1 int,
	window []int, overlap int) {
	
	if g0 == 0 && g1 == 0 {
		if yPtr != xPtr {
			copy(y[yPtr:yPtr+N], x[xPtr:xPtr+N])
		}
		return
	}
	
	gains := [3][3]int{
		{int(0.30664 * float32(1<<15)), int(0.21704 * float32(1<<15)), int(0.12964 * float32(1<<15))},
		{int(0.46387 * float32(1<<15)), int(0.26807 * float32(1<<15)), 0},
		{int(0.79980 * float32(1<<15)), int(0.10010 * float32(1<<15)), 0},
	}
	
	g00 := c.Mult16_16_P15(g0, gains[tapset0][0])
	g01 := c.Mult16_16_P15(g0, gains[tapset0][1])
	g02 := c.Mult16_16_P15(g0, gains[tapset0][2])
	g10 := c.Mult16_16_P15(g1, gains[tapset1][0])
	g11 := c.Mult16_16_P15(g1, gains[tapset1][1])
	g12 := c.Mult16_16_P15(g1, gains[tapset1][2])
	
	x1 := x[xPtr-T1+1]
	x2 := x[xPtr-T1]
	x3 := x[xPtr-T1-1]
	x4 := x[xPtr-T1-2]
	
	if g0 == g1 && T0 == T1 && tapset0 == tapset1 {
		overlap = 0
	}
	
	for i := 0; i < overlap; i++ {
		x0 := x[xPtr+i-T1+2]
		f := c.Mult16_16_Q15(window[i], window[i])
		invf := (1 << 15) - f
		
		term1 := c.Mult16_32_Q15(invf, g00) * x[xPtr+i-T0]
		term2 := c.Mult16_32_Q15(invf, g01) * (x[xPtr+i-T0+1] + x[xPtr+i-T0-1])
		term3 := c.Mult16_32_Q15(invf, g02) * (x[xPtr+i-T0+2] + x[xPtr+i-T0-2])
		term4 := c.Mult16_32_Q15(f, g10) * x2
		term5 := c.Mult16_32_Q15(f, g11) * (x1 + x3)
		term6 := c.Mult16_32_Q15(f, g12) * (x0 + x4)
		
		y[yPtr+i] = x[xPtr+i] + term1 + term2 + term3 + term4 + term5 + term6
		
		x4 = x3
		x3 = x2
		x2 = x1
		x1 = x0
	}
	
	if g1 == 0 {
		if yPtr != xPtr {
			copy(y[yPtr+overlap:yPtr+N], x[xPtr+overlap:xPtr+N])
		}
		return
	}
	
	c.CombFilterConst(y, yPtr+overlap, x, xPtr+overlap, T1, N-overlap, g10, g11, g12)
}

func (c *CeltCommon) InitCaps(mode *CeltMode, cap []int, LM int, channels int) {
	for i := 0; i < mode.nbEBands; i++ {
		N := (mode.eBands[i+1] - mode.eBands[i]) << LM
		idx := mode.nbEBands*(2*LM+channels-1) + i
		cap[i] = (mode.cache.caps[idx] + 64) * channels * N >> 2
	}
}

// Helper methods for fixed-point arithmetic
func (c *CeltCommon) Mult16_32_Q15(a, b int) int {
	return (a * b) >> 15
}

func (c *CeltCommon) Mult16_16_P15(a, b int) int {
	return (a * b) >> 15
}

func (c *CeltCommon) Mult16_16_Q15(a, b int) int {
	return (a * b) >> 15
}