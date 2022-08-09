package math

import "math"

type FloatingPoint interface {
	~float32 | ~float64
}

func abs[T FloatingPoint](in T) T {
	return T(math.Abs(float64(in)))
}
