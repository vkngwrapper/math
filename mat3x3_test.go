package math

import (
	"github.com/stretchr/testify/require"
	"math"
	"testing"
)

func TestMat3x3_SetMat3x3(t *testing.T) {
	var input Mat3x3[float32]
	axis := Vec3[float32]{X: 0.7071, Y: 0.7071, Z: 0}
	input.SetAxisAngle(&axis, math.Pi/2.0)

	var m Mat3x3[float32]
	m.SetMat3x3(&input)

	var outAxis Vec3[float32]
	var outAngle float64

	m.GetAxisAngle(&outAxis, &outAngle)

	require.InDelta(t, 0.7071, outAxis.X, 0.0001)
	require.InDelta(t, 0.7071, outAxis.Y, 0.0001)
	require.InDelta(t, 0, outAxis.Z, 0.0001)
	require.InDelta(t, math.Pi/2.0, outAngle, 0.0001)
}

func TestMat3x3_SetMat4x4(t *testing.T) {
	var input Mat4x4[float32]
	axis := Vec3[float32]{X: 0.7071, Y: 0.7071, Z: 0}
	input.SetAxisAngle(&axis, math.Pi/2.0)

	var m Mat3x3[float32]
	m.SetMat4x4(&input)

	var outAxis Vec3[float32]
	var outAngle float64

	m.GetAxisAngle(&outAxis, &outAngle)

	require.InDelta(t, 0.7071, outAxis.X, 0.0001)
	require.InDelta(t, 0.7071, outAxis.Y, 0.0001)
	require.InDelta(t, 0, outAxis.Z, 0.0001)
	require.InDelta(t, math.Pi/2.0, outAngle, 0.0001)
}

func TestMat3x3_SetMultMatrix3x3(t *testing.T) {
	var first Mat3x3[float32]
	first.SetRotationY(math.Pi / 2.0)

	var second Mat3x3[float32]
	second.SetRotationX(math.Pi / 2.0)

	var m3 Mat3x3[float32]
	m3.SetMultMatrix3x3(&second, &first)

	input := Vec3[float32]{X: 1, Y: 0, Z: 0}
	input.Transform(&m3)

	require.InDelta(t, 0, input.X, 0.0001)
	require.InDelta(t, 1, input.Y, 0.0001)
	require.InDelta(t, 0, input.Z, 0.0001)
}

func TestMat3x3_MultMatrix3x3(t *testing.T) {
	var first Mat3x3[float32]
	first.SetRotationY(math.Pi / 2.0)

	var m3 Mat3x3[float32]
	m3.SetRotationX(math.Pi / 2.0)

	m3.MultMatrix3x3(&first)

	input := Vec3[float32]{X: 1, Y: 0, Z: 0}
	input.Transform(&m3)

	require.InDelta(t, 0, input.X, 0.0001)
	require.InDelta(t, 1, input.Y, 0.0001)
	require.InDelta(t, 0, input.Z, 0.0001)
}

func TestMat3x3_SetApplyTransform(t *testing.T) {
	var first Mat3x3[float32]
	first.SetRotationY(math.Pi / 2.0)

	var second Mat3x3[float32]
	second.SetRotationX(math.Pi / 2.0)

	var m3 Mat3x3[float32]
	m3.SetApplyTransform(&first, &second)

	input := Vec3[float32]{X: 1, Y: 0, Z: 0}
	input.Transform(&m3)

	require.InDelta(t, 0, input.X, 0.0001)
	require.InDelta(t, 1, input.Y, 0.0001)
	require.InDelta(t, 0, input.Z, 0.0001)
}

func TestMat3x3_ApplyTransform(t *testing.T) {
	var m3 Mat3x3[float32]
	m3.SetRotationY(math.Pi / 2.0)

	var second Mat3x3[float32]
	second.SetRotationX(math.Pi / 2.0)
	m3.ApplyTransform(&second)

	input := Vec3[float32]{X: 1, Y: 0, Z: 0}
	input.Transform(&m3)

	require.InDelta(t, 0, input.X, 0.0001)
	require.InDelta(t, 1, input.Y, 0.0001)
	require.InDelta(t, 0, input.Z, 0.0001)
}

func TestMat3x3_SetRotationEulers(t *testing.T) {
	var eulers Mat3x3[float32]
	eulers.SetRotationEulers(math.Pi/2, -math.Pi/2, math.Pi)

	point := Vec3[float32]{-1, 0, 0}
	point.Transform(&eulers)

	require.InDelta(t, 0, point.X, 0.0001)
	require.InDelta(t, -1, point.Y, 0.0001)
	require.InDelta(t, 0, point.Z, 0.0001)
}

func TestMat3x3_RotateX(t *testing.T) {
	var rotation Mat3x3[float32]
	rotation.SetIdentity()
	rotation.RotateX(math.Pi / 2)

	point := Vec3[float32]{0, 1, 1}
	point.Transform(&rotation)

	require.InDelta(t, 0, point.X, 0.0001)
	require.InDelta(t, -1, point.Y, 0.0001)
	require.InDelta(t, 1, point.Z, 0.0001)
}

func TestMat3x3_RotateY(t *testing.T) {
	var rotation Mat3x3[float32]
	rotation.SetIdentity()
	rotation.RotateY(math.Pi / 2)

	point := Vec3[float32]{1, 0, 1}
	point.Transform(&rotation)

	require.InDelta(t, 1, point.X, 0.0001)
	require.InDelta(t, 0, point.Y, 0.0001)
	require.InDelta(t, -1, point.Z, 0.0001)
}

func TestMat3x3_RotateZ(t *testing.T) {
	var rotation Mat3x3[float32]
	rotation.SetIdentity()
	rotation.RotateZ(math.Pi / 2)

	point := Vec3[float32]{1, 1, 0}
	point.Transform(&rotation)

	require.InDelta(t, -1, point.X, 0.0001)
	require.InDelta(t, 1, point.Y, 0.0001)
	require.InDelta(t, 0, point.Z, 0.0001)
}

func TestMat3x3_SetRotationZ(t *testing.T) {
	var rotation Mat3x3[float32]
	rotation.SetRotationZ(math.Pi / 2)

	point := Vec3[float32]{1, 1, 0}
	point.Transform(&rotation)

	require.InDelta(t, -1, point.X, 0.0001)
	require.InDelta(t, 1, point.Y, 0.0001)
	require.InDelta(t, 0, point.Z, 0.0001)
}

func TestMat3x3_SetRotationAroundAxis(t *testing.T) {
	axis := Vec3[float32]{0, 1, 0}
	var rotation Mat3x3[float32]
	rotation.SetRotationAroundAxis(&axis, math.Pi)

	point := Vec3[float32]{1, 0, 1}
	point.Transform(&rotation)

	require.InDelta(t, -1, point.X, 0.0001)
	require.InDelta(t, 0, point.Y, 0.0001)
	require.InDelta(t, -1, point.Z, 0.0001)
}

func TestMat3x3_RotateAroundAxis(t *testing.T) {
	axis := Vec3[float32]{0, 1, 0}
	var rotation Mat3x3[float32]
	rotation.SetIdentity()
	rotation.RotateAroundAxis(&axis, math.Pi)

	point := Vec3[float32]{1, 0, 1}
	point.Transform(&rotation)

	require.InDelta(t, -1, point.X, 0.0001)
	require.InDelta(t, 0, point.Y, 0.0001)
	require.InDelta(t, -1, point.Z, 0.0001)
}

func TestMat3x3_Transpose(t *testing.T) {
	m := Mat3x3[float32]{
		{1, 2, 3},
		{4, 5, 6},
		{7, 8, 9},
	}

	m.Transpose()

	require.InDelta(t, 1, m[0][0], 0.0001)
	require.InDelta(t, 2, m[1][0], 0.0001)
	require.InDelta(t, 3, m[2][0], 0.0001)
	require.InDelta(t, 4, m[0][1], 0.0001)
	require.InDelta(t, 5, m[1][1], 0.0001)
	require.InDelta(t, 6, m[2][1], 0.0001)
	require.InDelta(t, 7, m[0][2], 0.0001)
	require.InDelta(t, 8, m[1][2], 0.0001)
	require.InDelta(t, 9, m[2][2], 0.0001)
}

func TestMat3x3_Inverse(t *testing.T) {
	mat1 := Mat3x3[float32]{
		{0.6, 0.2, 0.3},
		{0.2, 0.7, 0.5},
		{0.3, 0.5, 0.7},
	}
	var inverse Mat3x3[float32]
	inverse.SetMat3x3(&mat1).Inverse()
	var identity Mat3x3[float32]
	identity.SetIdentity()

	require.True(t, mat1.MultMatrix3x3(&inverse).Equal(&identity, 0.0001))
}

func TestMat3x3_Orthonormalize(t *testing.T) {
	m := Mat3x3[float32]{
		{3, 0.15, 0},
		{0, 1, 0},
		{0.05, 0, 1},
	}
	m.Orthonormalize()

	v1 := Vec3[float32]{X: m[0][0], Y: m[0][1], Z: m[0][2]}
	v2 := Vec3[float32]{X: m[1][0], Y: m[1][1], Z: m[1][2]}
	v3 := Vec3[float32]{X: m[2][0], Y: m[2][1], Z: m[2][2]}

	require.InDelta(t, 1, v1.LenSqr(), 0.0001)
	require.InDelta(t, 1, v2.LenSqr(), 0.0001)
	require.InDelta(t, 1, v3.LenSqr(), 0.0001)
	require.InDelta(t, 0, v1.DotProduct(&v2), 0.0001)
	require.InDelta(t, 0, v1.DotProduct(&v3), 0.0001)
	require.InDelta(t, 0, v2.DotProduct(&v3), 0.0001)
}

func TestMat3x3_SetShearX(t *testing.T) {
	var shear Mat3x3[float32]
	shear.SetShearX(0.5, 1)

	point := Vec3[float32]{1, 0, 0}
	point.Transform(&shear)

	require.InDelta(t, 1.0, point.X, 0.0001)
	require.InDelta(t, 0.5, point.Y, 0.0001)
	require.InDelta(t, 1.0, point.Z, 0.0001)
}

func TestMat3x3_ShearX_Identity(t *testing.T) {
	var shear Mat3x3[float32]
	shear.SetIdentity()
	shear.ShearX(0.5, 1)

	point := Vec3[float32]{1, 0, 0}
	point.Transform(&shear)

	require.InDelta(t, 1.0, point.X, 0.0001)
	require.InDelta(t, 0.5, point.Y, 0.0001)
	require.InDelta(t, 1.0, point.Z, 0.0001)
}

func TestMat3x3_SetShearY(t *testing.T) {
	var shear Mat3x3[float32]
	shear.SetShearY(0.5, 1)

	point := Vec3[float32]{0, 1, 0}
	point.Transform(&shear)

	require.InDelta(t, 0.5, point.X, 0.0001)
	require.InDelta(t, 1.0, point.Y, 0.0001)
	require.InDelta(t, 1.0, point.Z, 0.0001)
}

func TestMat3x3_ShearY_Identity(t *testing.T) {
	var shear Mat3x3[float32]
	shear.SetIdentity()
	shear.ShearY(0.5, 1)

	point := Vec3[float32]{0, 1, 0}
	point.Transform(&shear)

	require.InDelta(t, 0.5, point.X, 0.0001)
	require.InDelta(t, 1.0, point.Y, 0.0001)
	require.InDelta(t, 1.0, point.Z, 0.0001)
}

func TestMat3x3_SetShearZ(t *testing.T) {
	var shear Mat3x3[float32]
	shear.SetShearZ(0.5, 1)

	point := Vec3[float32]{0, 0, 1}
	point.Transform(&shear)

	require.InDelta(t, 0.5, point.X, 0.0001)
	require.InDelta(t, 1.0, point.Y, 0.0001)
	require.InDelta(t, 1.0, point.Z, 0.0001)
}

func TestMat3x3_ShearZ_Identity(t *testing.T) {
	var shear Mat3x3[float32]
	shear.SetIdentity()
	shear.ShearZ(0.5, 1)

	point := Vec3[float32]{0, 0, 1}
	point.Transform(&shear)

	require.InDelta(t, 0.5, point.X, 0.0001)
	require.InDelta(t, 1.0, point.Y, 0.0001)
	require.InDelta(t, 1.0, point.Z, 0.0001)
}

func TestMat3x3_Scale(t *testing.T) {
	var scale Mat3x3[float32]
	scale.SetIdentity()
	scale.Scale(1.5, 2, 3)

	point := Vec3[float32]{1, 2, 3}
	point.Transform(&scale)

	require.InDelta(t, 1.5, point.X, 0.0001)
	require.InDelta(t, 4, point.Y, 0.0001)
	require.InDelta(t, 9, point.Z, 0.0001)
}
