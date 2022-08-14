package math

import (
	"math"
)

//go:noescape
func __m4x4RotateY(m *Mat4x4, sin, cos float32)

func mat4x4RotateYSSE(m *Mat4x4, angleRad float32) {
	cos := float32(math.Cos(float64(angleRad)))
	sin := float32(math.Sin(float64(angleRad)))

	__m4x4RotateY(m, sin, cos)
}
