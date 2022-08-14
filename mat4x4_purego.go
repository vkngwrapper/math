package math

import "math"

func mat4x4RotateYPureGo(m *Mat4x4, angleRad float32) {
	cos := float32(math.Cos(float64(angleRad)))
	sin := float32(math.Sin(float64(angleRad)))

	m00 := m[0][0]*cos - m[2][0]*sin
	m01 := m[0][1]*cos - m[2][1]*sin
	m02 := m[0][2]*cos - m[2][2]*sin
	m03 := m[0][3]*cos - m[2][3]*sin

	m20 := m[0][0]*sin + m[2][0]*cos
	m21 := m[0][1]*sin + m[2][1]*cos
	m22 := m[0][2]*sin + m[2][2]*cos
	m23 := m[0][3]*sin + m[2][3]*cos

	m[0][0] = m00
	m[0][1] = m01
	m[0][2] = m02
	m[0][3] = m03
	m[2][0] = m20
	m[2][1] = m21
	m[2][2] = m22
	m[2][3] = m23
}
