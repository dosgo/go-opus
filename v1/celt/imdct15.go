package celt

import (
	"math"
)

const PI = math.Pi

type Complex complex64

func NewComplex(re, im float32) Complex {
	return Complex(complex(re, im))
}

func (c Complex) Real() float32 {
	return real(c)
}

func (c Complex) Imag() float32 {
	return imag(c)
}

func (c Complex) Conj() Complex {
	return Complex(complex(real(c), -imag(c)))
}

func (c Complex) Scale(s float32) Complex {
	return Complex(complex(real(c)*s, imag(c)*s))
}

type IMDCT15 struct {
	n       int
	len2    int
	len4    int
	tmp     []Complex
	exptab  [][]Complex
	twiddle []Complex
}

func p2len(p2 int) int {
	return 15 * (1 << uint(p2))
}

var FACT = []Complex{
	NewComplex(0.30901699437494745, 0.95105651629515353),
	NewComplex(-0.80901699437494734, 0.58778525229247325),
}

func mulc(a, b Complex) (float32, float32, float32, float32) {
	ar, ai := real(a), imag(a)
	br, bi := real(b), imag(b)
	return ar * br, ar * bi, ai * br, ai * bi
}

func m_c(inp Complex) [4]Complex {
	rr0, ri0, ir0, ii0 := mulc(inp, FACT[0])
	rr1, ri1, ir1, ii1 := mulc(inp, FACT[1])
	return [4]Complex{
		NewComplex(rr0-ii0, ir0+ri0),
		NewComplex(rr1-ii1, ir1+ri1),
		NewComplex(rr1+ii1, ir1-ri1),
		NewComplex(rr0+ii0, ir0-ri0),
	}
}

func fft5(inp []Complex, stride int) [5]Complex {
	z := [4][4]Complex{
		m_c(inp[1*stride]),
		m_c(inp[2*stride]),
		m_c(inp[3*stride]),
		m_c(inp[4*stride]),
	}

	return [5]Complex{
		inp[0] + inp[1*stride] + inp[2*stride] + inp[3*stride] + inp[4*stride],
		inp[0] + z[0][0] + z[1][1] + z[2][2] + z[3][3],
		inp[0] + z[0][1] + z[1][3] + z[2][0] + z[3][2],
		inp[0] + z[0][2] + z[1][0] + z[2][3] + z[3][1],
		inp[0] + z[0][3] + z[1][2] + z[2][1] + z[3][0],
	}
}

func NewIMDCT15(n int) *IMDCT15 {
	len2 := p2len(n)
	len := len2 * 2
	len4 := len2 / 2

	tmp := make([]Complex, len*2)

	// Create twiddle factors
	twiddle := make([]Complex, len2-len4)
	for i := range twiddle {
		idx := i + len4
		v := 2 * PI * (float32(idx) + 0.125) / float32(len)
		twiddle[i] = NewComplex(float32(math.Cos(float64(v))), float32(math.Sin(float64(v))))
	}

	// Create exponential tables
	exptab := make([][]Complex, 6)
	for i := range exptab {
		length := p2len(i)
		if length < 19 {
			length = 19
		}
		tab := make([]Complex, length)
		for j := range tab {
			v := 2 * PI * float32(j) / float32(length)
			tab[j] = NewComplex(float32(math.Cos(float64(v))), float32(math.Sin(float64(v))))
		}
		exptab[i] = tab
	}

	// Pad first exponential table
	for i := 0; i < 4; i++ {
		exptab[0] = append(exptab[0], exptab[0][i])
	}

	return &IMDCT15{
		n:       n,
		len2:    len2,
		len4:    len4,
		tmp:     tmp,
		exptab:  exptab,
		twiddle: twiddle,
	}
}

func (imdct *IMDCT15) fft15(out []Complex, inp []Complex, stride int) {
	exptab := imdct.exptab[0]

	tmp0 := fft5(inp, stride*3)
	tmp1 := fft5(inp[1*stride:], stride*3)
	tmp2 := fft5(inp[2*stride:], stride*3)

	for i := 0; i < 5; i++ {
		t0 := tmp0[i]
		t1 := tmp1[i]
		t2 := tmp2[i]

		e1 := t1 * exptab[i]
		e2 := t2 * exptab[2*i]
		out[i] = t0 + e1 + e2

		e1 = t1 * exptab[i+5]
		e2 = t2 * exptab[2*(i+5)]
		out[i+5] = t0 + e1 + e2

		e1 = t1 * exptab[i+10]
		e2 = t2 * exptab[2*i+5]
		out[i+10] = t0 + e1 + e2
	}
}

func (imdct *IMDCT15) fft_calc(n int, out []Complex, inp []Complex, stride int) {
	if n > 0 {
		exptab := imdct.exptab[n]
		len2 := p2len(n)

		imdct.fft_calc(n-1, out[:len2], inp, stride*2)
		imdct.fft_calc(n-1, out[len2:], inp[stride:], stride*2)

		for i := 0; i < len2; i++ {
			e := out[i+len2] * exptab[i]
			o := out[i]
			out[i] = o + e
			out[i+len2] = o - e
		}
	} else {
		imdct.fft15(out, inp, stride)
	}
}

func (imdct *IMDCT15) IMDCT15Half(out []float32, inp []float32, stride int, scale float32) {
	// Create a complex view of the output buffer
	dst := make([]Complex, len(out)/2)

	len8 := imdct.len4 / 2
	start := (imdct.len2 - 1) * stride

	// Fill tmp buffer with complex values
	for i := 0; i < imdct.len2; i++ {
		re := inp[start-2*stride*i]
		im := inp[2*stride*i]
		imdct.tmp[i] = NewComplex(re, im) * imdct.twiddle[i]
	}

	// Perform FFT calculation
	imdct.fft_calc(imdct.n, dst, imdct.tmp, 1)

	// Process the output
	for i := 0; i < len8; i++ {
		decr := len8 - i - 1
		incr := len8 + i

		re0im1 := NewComplex(dst[decr].Imag(), dst[decr].Real()) *
			NewComplex(imdct.twiddle[decr].Imag(), imdct.twiddle[decr].Imag())

		re1im0 := NewComplex(dst[incr].Imag(), dst[incr].Real()) *
			NewComplex(imdct.twiddle[incr].Imag(), imdct.twiddle[incr].Imag())

		dst[decr] = NewComplex(re0im1.Real(), re1im0.Imag()).Scale(scale)
		dst[incr] = NewComplex(re1im0.Real(), re0im1.Imag()).Scale(scale)
	}

	// Copy complex results back to real output
	for i, c := range dst {
		out[2*i] = real(c)
		out[2*i+1] = imag(c)
	}
}
