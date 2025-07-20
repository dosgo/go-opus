package celt

// mul16 乘法：int16 * int16 -> int32
func mul16(a, b int16) int32 {
	return int32(a) * int32(b)
}

// frac_mul16 定点小数乘法：(a * b) 右移15位 + 四舍五入
func frac_mul16(a, b int16) int16 {
	v := mul16(a, b)
	return int16((16384 + v) >> 15)
}

// cos 余弦近似计算（多项式逼近）
func cos(x int16) int16 {
	v := (mul16(x, x) + 4096) >> 13
	xv := int16(v)
	return 1 + (32767 - xv) + frac_mul16(xv, -7651+frac_mul16(xv, 8277+frac_mul16(-626, xv)))
}

// celtIlog2 计算整数的二进制对数（高位1的位置）
func celtIlog2(n int32) int {
	if n <= 0 {
		return 0
	}
	r := 0
	if n >= 1<<16 {
		n >>= 16
		r += 16
	}
	if n >= 256 {
		n >>= 8
		r += 8
	}
	if n >= 16 {
		n >>= 4
		r += 4
	}
	if n >= 4 {
		n >>= 2
		r += 2
	}
	if n >= 2 {
		r += 1
	}
	return r
}

// log2tan 计算 log2(tan) 的近似值
func log2tan(isin, icos int32) int32 {
	ls := celtIlog2(isin)
	lc := celtIlog2(icos)

	// 归一化到15位精度
	sin16 := int16(isin << (15 - ls))
	cos16 := int16(icos << (15 - lc))

	// 计算对数缩放值
	s := int32(ls<<11) + int32(frac_mul16(sin16, frac_mul16(sin16, -2597)+7932))
	c := int32(lc<<11) + int32(frac_mul16(cos16, frac_mul16(cos16, -2597)+7932))

	return s - c
}
