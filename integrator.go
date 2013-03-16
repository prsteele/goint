package goint

import (
	"math"
)

/* This package provides a one-dimentional numeric integrator based on
Newton-Cotes formulas. */

type Function func(x float64) float64

type Integrator func(f Function, a, b float64, nsteps int) float64

func SimpsonIntegration(f Function, a, b, err float64) float64 {
	// Get an initial estimate
	ret := (b - a) * (f(a) + 4*f((a+b)/2.0) + f(b)) / 6.0

	points := []float64{a, b}
	done := false
	for !done {
		// Get a refined estimate
		points = refinedPoints(points)

		refined := 0.0
		L := points[0]
		fL := f(L)
		for _, R := range points[1:] {
			M := (L + R) / 2.0
			fM := f(M)
			fR := f(R)

			refined += (R - L) * (fL + 4*fM + fR) / 6.0
			L = R
			fL = fR
		}

		if math.Abs(ret-refined) < err {
			done = true
		}

		ret = refined
	}

	return ret
}

func BooleIntegration(f Function, a, b, err float64) float64 {
	// Get an initial estimate
	ret := (b - a) * (f(a) + 4*f((a+b)/2.0) + f(b)) / 6.0

	points := []float64{a, b}
	done := false
	for !done {
		// Get a refined estimate
		points = refinedPoints(points)

		refined := 0.0
		L := points[0]
		fL := f(L)
		for _, R := range points[1:] {
			h := (R - L) / 4.0
			f2 := f(L + h)
			f3 := f(L + 2*h)
			f4 := f(L + 3*h)
			fR := f(R) // R == L + 4 * h

			refined += 2 * h * (7*fL + 32*f2 + 12*f3 + 32*f4 + 7*fR) / 45.0
			L = R
			fL = fR
		}

		if math.Abs(ret-refined) < err {
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
	refined := make([]float64, len(points)*2-1)

	ndx := 0
	for i := 0; i < len(points)-1; i++ {
		refined[ndx] = points[i]
		refined[ndx+1] = (points[i] + points[i+1]) / 2.0
		ndx += 2
	}
	refined[ndx] = points[len(points)-1]

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
		for x < b {
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
