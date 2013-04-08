package goint

import (
	"math"
)

/* This package provides a one-dimentional numeric integrator based on
Newton-Cotes formulas. */

type Function func(x float64) float64

type Integrator func(f Function, a, b, err float64) float64

func boolesrule(f Function, a, b float64) float64 {
	h := (b - a) / 4.0
	fa := f(a)
	f2 := f(a + h)
	f3 := f(a + 2*h)
	f4 := f(a + 3*h)
	fb := f(b) // b == L + 4 * h

	return 2 * h * (7*fa + 32*f2 + 12*f3 + 32*f4 + 7*fb) / 45.0
}

/* Integrate a function f over the interval [a, b] to within err. Both
/* a and b can be infinite. Integration will be done using Boole's
/* rule. */
func Integrate(f Function, a, b, err float64) float64 {
	var ret float64

	// Get an initial estimate, being conservative when there are infinities
	if math.IsInf(a, -1) || math.IsInf(b, 1) {
		ret = math.Inf(1)
	} else {
		ret = boolesrule(f, a, b)
	}

	points := []float64{a, b}
	done := false
	for !done {
		// Get a refined estimate
		points = refinedPoints(points)

		// Skip extreme points
		start := 1
		end := len(points)

		if math.IsInf(points[0], -1) {
			start += 1
		}

		if math.IsInf(points[end-1], 1) {
			end -= 1
		}

		refined := 0.0
		L := points[start-1]
		for _, R := range points[start:end] {
			refined += boolesrule(f, L, R)
			L = R
		}

		// Check for unbounded integrals
		if math.IsInf(ret, 1) && math.IsInf(refined, 1) {
			return ret
		} else if math.IsInf(ret, -1) && math.IsInf(refined, -1) {
			return ret
		} else if math.Abs(ret-refined) < err {
			done = true
		}

		ret = refined
	}

	return ret
}

/* Returns a new slice of values containing all the values in points
/* as well as the midpoint of each sequential pair in points. For example,
/*
/*   refinedPoints([]float64{0, 2, 4}) == []float64{0, 1, 2, 3, 4} */
func refinedPoints(points []float64) []float64 {
	// Check for infinite extremes with only two points specified
	if len(points) == 2 {
		if math.IsInf(points[0], -1) && math.IsInf(points[1], 1) {
			return []float64{points[0], 0, points[1]}
		} else if math.IsInf(points[0], -1) && points[1] >= 0 {
			return []float64{points[0], -1, points[1]}
		} else if math.IsInf(points[1], 1) && points[0] <= 0 {
			return []float64{points[0], 1, points[1]}
		}
	}

	refined := make([]float64, len(points)*2-1)

	// Check the left endpoint for -Inf
	if math.IsInf(points[0], -1) {
		refined[0] = points[0]
		refined[1] = points[1] * 2
	} else {
		refined[0] = points[0]
		refined[1] = (points[0] + points[1]) / 2
	}

	points_end := len(points) - 1
	refined_end := len(refined) - 1

	// Interpolate all intermediate points
	ndx := 2
	for i := 1; i < points_end-1; i++ {
		refined[ndx] = points[i]
		refined[ndx+1] = (points[i] + points[i+1]) / 2
		ndx += 2
	}

	// Check the right endpoint for +Inf
	if math.IsInf(points[points_end], 1) {
		refined[refined_end] = points[points_end]
		refined[refined_end-1] = points[points_end-1] * 2
		refined[refined_end-2] = points[points_end-1]
	} else {
		refined[refined_end] = points[points_end]
		refined[refined_end-1] = (points[points_end] + points[points_end-1]) / 2
		refined[refined_end-2] = points[points_end-1]
	}

	return refined
}
