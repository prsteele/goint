package main

import (
	"fmt"
	"goint"
	"math"
)

func f(x float64) float64 {
	// x^5 - x^4 + x^2 - 1
	xsq := x * x
	return xsq*(xsq*x-xsq+1) - 1
}

func main() {
	fmt.Println(goint.SimpsonIntegration(f, 0, 5, 1e-4))
	fmt.Println(goint.SimpsonIntegration(math.Exp, math.Inf(-1), 0, 1e-6))
}
