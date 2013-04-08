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

/* Sends float64 values along the provided channel spaced across the
/* provided interval. In all cases the values are in increasing order.
/*
/* If both a and b are finite, this function sends points along the
/* interval [a, b] with h between each point; the final interval may
/* be slightly smaller, as we always include the enpoints a and b. 
/*
/* If a == -Inf and b is finite, we compute a series of points 0, -h,
/* -h*s, -h*s^2, ... -h*s^q such that b - h * s^q < -1e50, where
/* s=1.1. We then return these points b - h * s^q, b - h * s^(q - 1),
/* ..., b along the channel.
/*
/* If a is finite and b == Inf, we compute a series of points 0, h,
/* h*s, h*s^2, ..., h*s^q such that a+h*s^q > 1e50 where s=1.1. We
/* then send the points a, a+h, a+h*s, a+h*s^2, ..., a+h*s^q along the
/* channel.
/*
/* Finally, if a == -Inf and b == Inf, we first send the points that
/* would be sent from the interval (a, 0], and then those from (0,
/* b). 
/*
/* Since we use an exponential spacing scheme for infinite endpoints,
/* more accurate results can be obtained by summing two (or more)
/* separate calls to the points function, where regions of interest
/* (e.g., regions with large derivatives) are computed over finite
/* regions, and unintersting (e.g., well-behaved) tails are called
/* over infinite regions. For example, consider integrating e^-x over
/* [0, Inf). We might instead sum the integrals over [0, 10] and [10,
/* Inf) to increase the number of points used in the region with 1 -
/* e^(-10) of the total mass. */
func points(a, b, h float64, c chan float64) {
	defer close(c)
	const Large = 1e50
	const Scale = 1.1

	if !math.IsInf(a, -1) && !math.IsInf(b, 1) {
		// Both endpoints are finite
		x := a
		for x <= b {
			c <- x
			x += h
		}
	} else if math.IsInf(1, -1) && math.IsInf(b, 1) {
		// Both endpoints are infinite

		// Find the smallest value we send
		x := -h
		for x > -Large {
			x *= Scale
		}

		// Send values from the smallest up to zero
		for x < -h {
			c <- x
			x /= Scale
		}
		c <- 0.0

		// Send values from h out to the largest values
		x = h
		for x < Large {
			c <- x
			x *= Scale
		}
	} else if math.IsInf(1, -1) {
		// Only the left endpoint is infinite. Find the smallest value we
		// send; in this case, x is an offset from b
		x := -h
		for b+x > -Large {
			x *= Scale
		}

		// Send values from the smallest up to zero
		for x < -h {
			c <- x + b
			x /= Scale
		}
		c <- b
	} else {
		// Only the left endpoint is infinite. In this case, x is an
		// offset from a
		x := h
		c <- a
		for a+x < Large {
			c <- a + x
			x *= Scale
		}
	}
}
