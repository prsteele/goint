package goint

import (
	"fmt"
	"math"
	"testing"
)

/* Computes the integral of f over [a, b], and determines if the
/* result is within h of the correct value. Returns (msg, ok), where
/* ok is true if there was an acceptable amount of error. If ok is
/* false, msg contains a string describing the error, otherwise it is
/* empty. */
func test_integral(f Function, a, b, h, correct float64) (string, bool) {
	int := Integrate(f, a, b, h)
	err := math.Abs(int - correct)

	if err > h {
		msg := fmt.Sprintf("%.3g differs from %.3g by more than %.3g", int, correct, h)
		return msg, false
	}

	return "", true
}

/* Test integrating e^x over infinite domains */
func TestExponential(t *testing.T) {
	const (
		h = 1e-8
	)

	f := math.Exp

	// Check (-Inf, 0]; should be 1
	if msg, ok := test_integral(f, math.Inf(-1), 0, h, 1); !ok {
		t.Error(msg)
	}

	// Check [0, +Inf); should be +Inf
	if msg, ok := test_integral(f, 0, math.Inf(1), h, math.Inf(1)); !ok {
		t.Error(msg)
	}

	// Now check that these results hold for -f
	f = func(x float64) float64 { return -math.Exp(x) }

	// Check (-Inf, 0]; should be -1
	if msg, ok := test_integral(f, math.Inf(-1), 0, h, -1); !ok {
		t.Error(msg)
	}

	// Check [0, +Inf); should be -Inf
	if msg, ok := test_integral(f, 0, math.Inf(1), h, math.Inf(-1)); !ok {
		t.Error(msg)
	}
}

/* Test integrating e^(-x) over infinite domains*/
func TestNegativeExponential(t *testing.T) {
	const (
		h       = 1e-8
		correct = 1
	)

	f := func(x float64) float64 { return math.Exp(-x) }

	// Check (-Inf, 0]; should be +Inf
	if msg, ok := test_integral(f, math.Inf(-1), 0, h, math.Inf(1)); !ok {
		t.Error(msg)
	}

	// Check [0, +Inf); should be 1
	if msg, ok := test_integral(f, 0, math.Inf(1), h, 1); !ok {
		t.Error(msg)
	}

	// Now check that these results hold for -f
	f = func(x float64) float64 { return -math.Exp(-x) }

	// Check (-Inf, 0]; should be -Inf
	if msg, ok := test_integral(f, math.Inf(-1), 0, h, math.Inf(-1)); !ok {
		t.Error(msg)
	}

	// Check [0, +Inf); should be -1
	if msg, ok := test_integral(f, 0, math.Inf(1), h, -1); !ok {
		t.Error(msg)
	}
}

func TestNormal(t *testing.T) {
	f := func(x float64) float64 {
		return math.Exp(-x*x/2) / math.Sqrt(2*math.Pi)
	}

	if msg, ok := test_integral(f, math.Inf(-1), 0, 1e-6, .5); !ok {
		t.Error(msg)
	}
}
