package math

import "math"

type FloatingPoint interface {
	~float32 | ~float64
}

func abs[T FloatingPoint](in T) T {
	return T(math.Abs(float64(in)))
}

func cosOneOverTwo[T FloatingPoint]() T {
	return 0.877582561890372716130286068203503191
}

func clamp[T FloatingPoint](v, min, max T) T {
	if v < min {
		return min
	}
	if v > max {
		return max
	}

	return v
}

func mix[T FloatingPoint](x, y, delta T) T {
	return x*(1.0-delta) + y*delta
}
