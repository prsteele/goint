package goint

import (
	"math"
	"testing"
)

/* Return a function representing the polynomial in x given by
/* coefs[0] + coefs[1] * x + ... */
func makePolynomial(coefs []float64) Function {
	return func(x float64) float64 {
		ret := 0.0
		xx := 1.0
		for _, c := range coefs {
			ret += xx * c
			xx *= x
		}

		return ret
	}
}

func polynomials() ([]Function, []Function) {
	// Polynomials to test
	Ps := []Function{
		makePolynomial([]float64{1}),                         // 1
		makePolynomial([]float64{1, 1}),                      // 1 + x
		makePolynomial([]float64{1, -1, 3}),                  // 1 - x + 3 x^2
		makePolynomial([]float64{0, 0, 0, 4}),                // 4 x^3
		makePolynomial([]float64{1, 0, 0, 0, -5}),            // 1 - 5 x^4
		makePolynomial([]float64{0, 0, 0, 0, 0, 0, 0, 0, 9}), // 9 x^8
	}

	// Their exact integrals (up to a constant)
	Is := []Function{
		makePolynomial([]float64{0, 1}),                         // x
		makePolynomial([]float64{0, 1, .5}),                     // x + .5 x^2
		makePolynomial([]float64{0, 1, -.5, 1}),                 // x - .5 x^2 + x^3
		makePolynomial([]float64{0, 0, 0, 0, 1}),                // x^4
		makePolynomial([]float64{0, 1, 0, 0, 0, -1}),            // x - x^5
		makePolynomial([]float64{0, 0, 0, 0, 0, 0, 0, 0, 0, 1}), // x^9
	}

	return Ps, Is
}

func TestPolynomials(t *testing.T) {
	const (
		h   = 1e-5
		a   = -1
		b   = 3
		err = 1e-7
	)

	Ps, Is := polynomials()

	for i, _ := range Ps {
		p := Ps[i]
		p_int := Is[i]

		computed_val := Integrate(p, a, b, h)
		correct_val := p_int(b) - p_int(a)
		computed_err := math.Abs(computed_val - correct_val)

		if computed_err > err {
			t.Errorf("Error %.3g exceeds acceptable error %.3g", computed_err, err)
		}
	}
}
