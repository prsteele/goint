package goint

import (
	"errors"
	"math"
)

/* This package provides a one-dimentional numeric integrator based on
Newton-Cotes formulas. */

type Function func(x float64) float64

type Integrator func(f Function, a, b float64, nsteps int) float64

func SimpsonIntegration(f Function, a, b, h float64) float64 {
	if h <= 0 {
		err := errors.New("Must specify a step width")
		panic(err)
	}

	// Get a channel of points ready
	pts := make(chan float64)
	go points(a, b, h, pts)

	// Compute Simpson's rule over each interval
	result := 0.0
	L := <-pts
	fL := f(L)
	for R := range pts {
		coef := (R - L) / 6.0
		fM := f((L + R) / 2.0)
		fR := f(R)

		result += coef * (fL + 4*fM + fR)
		L = R
		fL = fR
	}

	return result
}

func BooleIntegration(f Function, a, b, h float64) float64 {
	if h <= 0 {
		err := errors.New("Must specify a step width")
		panic(err)
	}

	// Get a channel of points ready
	pts := make(chan float64)
	go points(a, b, h, pts)

	// Compute Boole's rule over each interval
	result := 0.0
	L := <-pts
	fL := f(L)

	for R := range pts {
		w := (R - L) / 4.0
		coef := 2 * w / 45.0

		f2 := f(L + w)
		f3 := f(L + 2*w)
		f4 := f(L + 3*w)
		fR := f(R) // R == L + 4 * w

		result += coef * (7*fL + 32*f2 + 12*f3 + 32*f4 + 7*fR)
		L = R
		fL = fR
	}

	return result
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
