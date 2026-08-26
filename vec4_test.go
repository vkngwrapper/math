package math

import (
	"math"
	"testing"

	"github.com/stretchr/testify/require"
)

func TestVec4_SetVec2(t *testing.T) {
	v2 := Vec2{1, 3}

	var v4 Vec4
	v4.SetVec2(&v2)

	require.InDelta(t, 1, v4.X(), 0.0001)
	require.InDelta(t, 3, v4.Y(), 0.0001)
	require.InDelta(t, 0, v4.Z(), 0.0001)
	require.InDelta(t, 1, v4.W(), 0.0001)
}

func TestVec4_SetVec3(t *testing.T) {
	v3 := Vec3{1, 3, 5}

	var v4 Vec4
	v4.SetVec3(&v3)

	require.InDelta(t, 1, v4.X(), 0.0001)
	require.InDelta(t, 3, v4.Y(), 0.0001)
	require.InDelta(t, 5, v4.Z(), 0.0001)
	require.InDelta(t, 1, v4.W(), 0.0001)
}

func TestVec4_SetVec4(t *testing.T) {
	input := Vec4{1, 3, 5, 7}

	var v4 Vec4
	v4.SetVec4(&input)

	require.InDelta(t, 1, input.X(), 0.0001)
	require.InDelta(t, 3, input.Y(), 0.0001)
	require.InDelta(t, 5, input.Z(), 0.0001)
	require.InDelta(t, 7, input.W(), 0.0001)
}

func TestVec4_SetTransform(t *testing.T) {
	var mat Mat4x4
	mat.SetScale(2, 2, 2)
	mat.Translate(1, 1, 1)
	input := Vec4{1, 3, 5, 1}

	var v4 Vec4
	v4.SetTransform(&input, &mat)

	require.InDelta(t, 3, v4.X(), 0.0001)
	require.InDelta(t, 7, v4.Y(), 0.0001)
	require.InDelta(t, 11, v4.Z(), 0.0001)
	require.InDelta(t, 1, v4.W(), 0.0001)
}

func TestVec4_SetTransformHomogenous(t *testing.T) {
	var mat Mat4x4
	mat.SetScale(2, 2, 2)
	input := Vec4{1, 3, 5, 2}

	var v4 Vec4
	v4.SetTransformHomogenous(&input, &mat)

	require.InDelta(t, 1, v4.X(), 0.0001)
	require.InDelta(t, 3, v4.Y(), 0.0001)
	require.InDelta(t, 5, v4.Z(), 0.0001)
	require.InDelta(t, 1, v4.W(), 0.0001)
}

func TestVec4_SetRotateWithQuaternion(t *testing.T) {
	var quaternion Quaternion
	quaternion.SetRotationY(math.Pi)

	input := Vec4{1, 3, 5, 1}

	var v4 Vec4
	v4.SetRotateWithQuaternion(&input, &quaternion)

	require.InDelta(t, -1, v4.X(), 0.0001)
	require.InDelta(t, 3, v4.Y(), 0.0001)
	require.InDelta(t, -5, v4.Z(), 0.0001)
	require.InDelta(t, 1, v4.W(), 0.0001)
}

func TestVec4_SetRotate(t *testing.T) {
	normal := Vec3{0, 1, 0}
	input := Vec4{1, 3, 5, 1}

	var v4 Vec4
	v4.SetRotate(&input, math.Pi, &normal)

	require.InDelta(t, -1, v4.X(), 0.0001)
	require.InDelta(t, 3, v4.Y(), 0.0001)
	require.InDelta(t, -5, v4.Z(), 0.0001)
	require.InDelta(t, 1, v4.W(), 0.0001)
}

func TestVec4_Rotate(t *testing.T) {
	normal := Vec3{0, 1, 0}
	v4 := Vec4{1, 3, 5, 1}
	v4.Rotate(math.Pi, &normal)

	require.InDelta(t, -1, v4.X(), 0.0001)
	require.InDelta(t, 3, v4.Y(), 0.0001)
	require.InDelta(t, -5, v4.Z(), 0.0001)
	require.InDelta(t, 1, v4.W(), 0.0001)
}

func TestVec4_SetRotateX(t *testing.T) {
	input := Vec4{1, 3, 5, 1}

	var v4 Vec4
	v4.SetRotateX(&input, math.Pi)

	require.InDelta(t, 1, v4.X(), 0.0001)
	require.InDelta(t, -3, v4.Y(), 0.0001)
	require.InDelta(t, -5, v4.Z(), 0.0001)
	require.InDelta(t, 1, v4.W(), 0.0001)
}

func TestVec4_RotateX(t *testing.T) {
	v4 := Vec4{1, 3, 5, 1}
	v4.RotateX(math.Pi)

	require.InDelta(t, 1, v4.X(), 0.0001)
	require.InDelta(t, -3, v4.Y(), 0.0001)
	require.InDelta(t, -5, v4.Z(), 0.0001)
	require.InDelta(t, 1, v4.W(), 0.0001)
}

func TestVec4_SetRotateY(t *testing.T) {
	input := Vec4{1, 3, 5, 1}

	var v4 Vec4
	v4.SetRotateY(&input, math.Pi)

	require.InDelta(t, -1, v4.X(), 0.0001)
	require.InDelta(t, 3, v4.Y(), 0.0001)
	require.InDelta(t, -5, v4.Z(), 0.0001)
	require.InDelta(t, 1, v4.W(), 0.0001)
}

func TestVec4_RotateY(t *testing.T) {
	v4 := Vec4{1, 3, 5, 1}
	v4.RotateY(math.Pi)

	require.InDelta(t, -1, v4.X(), 0.0001)
	require.InDelta(t, 3, v4.Y(), 0.0001)
	require.InDelta(t, -5, v4.Z(), 0.0001)
	require.InDelta(t, 1, v4.W(), 0.0001)
}

func TestVec4_SetRotateZ(t *testing.T) {
	input := Vec4{1, 3, 5, 1}

	var v4 Vec4
	v4.SetRotateZ(&input, math.Pi)

	require.InDelta(t, -1, v4.X(), 0.0001)
	require.InDelta(t, -3, v4.Y(), 0.0001)
	require.InDelta(t, 5, v4.Z(), 0.0001)
	require.InDelta(t, 1, v4.W(), 0.0001)
}

func TestVec4_RotateZ(t *testing.T) {
	v4 := Vec4{1, 3, 5, 1}
	v4.RotateZ(math.Pi)

	require.InDelta(t, -1, v4.X(), 0.0001)
	require.InDelta(t, -3, v4.Y(), 0.0001)
	require.InDelta(t, 5, v4.Z(), 0.0001)
	require.InDelta(t, 1, v4.W(), 0.0001)
}

func TestVec4_SetNormalizeVec4(t *testing.T) {
	input := Vec4{1, 3, 5, 1}

	var v4 Vec4
	v4.SetNormalizeVec4(&input)

	require.InDelta(t, 1.0/6.0, v4.X(), 0.0001)
	require.InDelta(t, 0.5, v4.Y(), 0.0001)
	require.InDelta(t, 5.0/6.0, v4.Z(), 0.0001)
	require.InDelta(t, 1.0/6.0, v4.W(), 0.0001)
}

func TestVec4_SetAddVec4(t *testing.T) {
	a := Vec4{1, 3, 5, 1}
	b := Vec4{2, 4, 6, 2}

	var v4 Vec4
	v4.SetAddVec4(&a, &b)

	require.InDelta(t, 3, v4.X(), 0.0001)
	require.InDelta(t, 7, v4.Y(), 0.0001)
	require.InDelta(t, 11, v4.Z(), 0.0001)
	require.InDelta(t, 3, v4.W(), 0.0001)
}

func TestVec4_SetSubtractVec4(t *testing.T) {
	a := Vec4{1, 3, 5, 1}
	b := Vec4{2, 4, 6, 2}

	var v4 Vec4
	v4.SetSubtractVec4(&b, &a)

	require.InDelta(t, 1, v4.X(), 0.0001)
	require.InDelta(t, 1, v4.Y(), 0.0001)
	require.InDelta(t, 1, v4.Z(), 0.0001)
	require.InDelta(t, 1, v4.W(), 0.0001)
}

func TestVec4_SetScaleVec4(t *testing.T) {
	input := Vec4{1, 3, 5, 1}

	var v4 Vec4
	v4.SetScale(&input, 2.5)

	require.InDelta(t, 2.5, v4.X(), 0.0001)
	require.InDelta(t, 7.5, v4.Y(), 0.0001)
	require.InDelta(t, 12.5, v4.Z(), 0.0001)
	require.InDelta(t, 2.5, v4.W(), 0.0001)
}

func TestVec4_SetLerpVec4(t *testing.T) {
	a := Vec4{1, 3, 5, 1}
	b := Vec4{2, 4, 6, 2}

	var v4 Vec4
	v4.SetLerp(&a, &b, 0.5)

	require.InDelta(t, 1.5, v4.X(), 0.0001)
	require.InDelta(t, 3.5, v4.Y(), 0.0001)
	require.InDelta(t, 5.5, v4.Z(), 0.0001)
	require.InDelta(t, 1.5, v4.W(), 0.0001)
}

func TestVec4_LerpVec4(t *testing.T) {
	left := Vec4{1, 3, 5, 1}
	right := Vec4{2, 4, 6, 2}

	left.Lerp(&right, 0.5)

	require.InDelta(t, 1.5, left.X(), 0.0001)
	require.InDelta(t, 3.5, left.Y(), 0.0001)
	require.InDelta(t, 5.5, left.Z(), 0.0001)
	require.InDelta(t, 1.5, left.W(), 0.0001)
}

func TestVec4_Normalize(t *testing.T) {
	v4 := Vec4{1, 3, 5, 1}
	v4.Normalize()

	require.InDelta(t, 1.0/6.0, v4.X(), 0.0001)
	require.InDelta(t, 0.5, v4.Y(), 0.0001)
	require.InDelta(t, 5.0/6.0, v4.Z(), 0.0001)
	require.InDelta(t, 1.0/6.0, v4.W(), 0.0001)
}

func TestVec4_LenSqr(t *testing.T) {
	v4 := Vec4{1, 3, 5, 1}
	require.InDelta(t, 36, v4.LenSqr(), 0.0001)
}

func TestVec4_AddVec4(t *testing.T) {
	left := Vec4{1, 3, 5, 7}
	right := Vec4{2, 4, 6, 8}

	left.AddVec4(&right)

	require.InDelta(t, 3, left.X(), 0.0001)
	require.InDelta(t, 7, left.Y(), 0.0001)
	require.InDelta(t, 11, left.Z(), 0.0001)
	require.InDelta(t, 15, left.W(), 0.0001)
}

func TestVec4_SubtractVec4(t *testing.T) {
	left := Vec4{2, 4, 6, 8}
	right := Vec4{1, 3, 5, 7}

	left.SubtractVec4(&right)

	require.InDelta(t, 1, left.X(), 0.0001)
	require.InDelta(t, 1, left.Y(), 0.0001)
	require.InDelta(t, 1, left.Z(), 0.0001)
	require.InDelta(t, 1, left.W(), 0.0001)
}

func TestVec4_DotProduct(t *testing.T) {
	left := Vec4{1, 3, 5, 7}
	right := Vec4{2, 4, 6, 8}

	require.InDelta(t, 100, left.DotProduct(&right), 0.0001)
}

func TestVec4_Scale(t *testing.T) {
	v4 := Vec4{1, 3, 5, 7}
	v4.Scale(2)

	require.InDelta(t, 2, v4.X(), 0.0001)
	require.InDelta(t, 6, v4.Y(), 0.0001)
	require.InDelta(t, 10, v4.Z(), 0.0001)
	require.InDelta(t, 14, v4.W(), 0.0001)
}

func TestVec4_TransformHomogenous(t *testing.T) {
	var mat Mat4x4
	mat.SetScale(2, 2, 2)
	v4 := Vec4{1, 3, 5, 2}
	v4.TransformHomogenous(&mat)

	require.InDelta(t, 1, v4.X(), 0.0001)
	require.InDelta(t, 3, v4.Y(), 0.0001)
	require.InDelta(t, 5, v4.Z(), 0.0001)
	require.InDelta(t, 1, v4.W(), 0.0001)
}
