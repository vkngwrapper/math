//go:build !amd64

package math

func (m *Mat4x4) RotateY(angleRad float32) *Mat4x4 {
	mat4x4RotateYPureGo(m, angleRad)
	return m
}
