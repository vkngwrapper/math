package math

import "math"

func abs(in float32) float32 {
	return float32(math.Abs(float64(in)))
}

func cosOneOverTwo() float32 {
	return 0.877582561890372716130286068203503191
}

func clamp(v, min, max float32) float32 {
	if v < min {
		return min
	}
	if v > max {
		return max
	}

	return v
}

func mix(x, y, delta float32) float32 {
	return x*(1.0-delta) + y*delta
}
