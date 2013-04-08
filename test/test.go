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
	fmt.Println(goint.Integrate(f, 0, 5, 1e-4))
	fmt.Println(goint.Integrate(math.Exp, math.Inf(-1), 0, 1e-6))
	fmt.Println(goint.Integrate(func(x float64) float64 { return math.Exp(-x) }, 0, math.Inf(1), 1e-6))
}
