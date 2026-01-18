package math

import "math"

func ToRadians[T FloatingPoint](degrees T) float64 {
	return float64(degrees) * math.Pi / 180.0
}

func ToDegrees[T FloatingPoint](radians float64) T {
	return T(radians * 180.0 / math.Pi)
}
