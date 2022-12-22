package math

import (
	"github.com/stretchr/testify/require"
	"math"
	"testing"
)

func TestVec3_SetVec2(t *testing.T) {
	v2 := Vec2[float32]{X: 3, Y: 7}

	var vec Vec3[float32]
	vec.SetVec2(&v2)

	require.InDelta(t, 3.0, vec.X, 0.0001)
	require.InDelta(t, 7.0, vec.Y, 0.0001)
	require.InDelta(t, 0.0, vec.Z, 0.0001)
}

func TestVec3_SetVec4(t *testing.T) {
	v4 := Vec4[float32]{X: 3, Y: 7, Z: 11, W: 13}

	var vec Vec3[float32]
	vec.SetVec4(&v4)

	require.InDelta(t, 3.0, vec.X, 0.0001)
	require.InDelta(t, 7.0, vec.Y, 0.0001)
	require.InDelta(t, 11.0, vec.Z, 0.0001)
}

func TestVec3_SetHomogenousVec4(t *testing.T) {
	v4 := Vec4[float32]{X: 3, Y: 7, Z: 11, W: 2}

	var vec Vec3[float32]
	vec.SetHomogenousVec4(&v4)

	require.InDelta(t, 1.5, vec.X, 0.0001)
	require.InDelta(t, 3.5, vec.Y, 0.0001)
	require.InDelta(t, 5.5, vec.Z, 0.0001)
}

func TestVec3_SubtractVec3(t *testing.T) {
	v1 := Vec3[float32]{X: 5, Y: 9, Z: 13}
	v2 := Vec3[float32]{X: 1, Y: 2, Z: 4}

	var vec Vec3[float32]
	vec.SetSubtractVec3(&v1, &v2)

	require.InDelta(t, 4.0, vec.X, 0.0001)
	require.InDelta(t, 7.0, vec.Y, 0.0001)
	require.InDelta(t, 9.0, vec.Z, 0.0001)
}

func TestVec3_SetAddVec3(t *testing.T) {
	v1 := Vec3[float32]{X: 5, Y: 9, Z: 13}
	v2 := Vec3[float32]{X: 1, Y: 2, Z: 4}

	var vec Vec3[float32]
	vec.SetAddVec3(&v1, &v2)

	require.InDelta(t, 6.0, vec.X, 0.0001)
	require.InDelta(t, 11.0, vec.Y, 0.0001)
	require.InDelta(t, 17.0, vec.Z, 0.0001)
}

func TestVec3_SetTransform(t *testing.T) {
	v1 := Vec3[float32]{X: 1, Y: 2, Z: 3}
	var mat Mat3x3[float32]
	mat.SetScale(2, 2, 2)

	var vec Vec3[float32]
	vec.SetTransform(&v1, &mat)

	require.InDelta(t, 2.0, vec.X, 0.0001)
	require.InDelta(t, 4.0, vec.Y, 0.0001)
	require.InDelta(t, 6.0, vec.Z, 0.0001)
}

func TestVec3_SetTransformHomogenous(t *testing.T) {
	v1 := Vec3[float32]{X: 1, Y: 2, Z: 3}
	var mat Mat4x4[float32]
	mat.SetScale(2, 2, 2)
	mat.Translate(1.5, 3.4, 5.3)
	mat.Scale(2, 2, 2)

	var vec Vec3[float32]
	vec.SetTransformHomogenous(&v1, &mat)

	require.InDelta(t, 7, vec.X, 0.0001)
	require.InDelta(t, 14.8, vec.Y, 0.0001)
	require.InDelta(t, 22.6, vec.Z, 0.0001)
}

func TestVec3_TransformHomogenous(t *testing.T) {
	v1 := Vec3[float32]{X: 1, Y: 2, Z: 3}
	var mat Mat4x4[float32]
	mat.SetScale(2, 2, 2)
	mat.Translate(1.5, 3.4, 5.3)
	mat.Scale(2, 2, 2)

	v1.TransformHomogenous(&mat)

	require.InDelta(t, 7, v1.X, 0.0001)
	require.InDelta(t, 14.8, v1.Y, 0.0001)
	require.InDelta(t, 22.6, v1.Z, 0.0001)
}

func TestVec3_SetRotate(t *testing.T) {
	axis := Vec3[float32]{X: 0, Y: 1, Z: 0}
	input := Vec3[float32]{X: 1, Y: 0, Z: 0}

	var vec Vec3[float32]
	vec.SetRotate(&input, math.Pi, &axis)

	require.InDelta(t, -1.0, vec.X, 0.0001)
	require.InDelta(t, 0, vec.Y, 0.0001)
	require.InDelta(t, 0, vec.Z, 0.0001)
}

func TestVec3_SetRotateX(t *testing.T) {
	input := Vec3[float32]{X: 0, Y: 1, Z: 0}

	var vec Vec3[float32]
	vec.SetRotateX(&input, math.Pi/2.0)

	require.InDelta(t, 0, vec.X, 0.0001)
	require.InDelta(t, 0, vec.Y, 0.0001)
	require.InDelta(t, 1, vec.Z, 0.0001)
}

func TestVec3_SetRotateY(t *testing.T) {
	input := Vec3[float32]{X: 1, Y: 0, Z: 0}

	var vec Vec3[float32]
	vec.SetRotateY(&input, math.Pi/2.0)

	require.InDelta(t, 0, vec.X, 0.0001)
	require.InDelta(t, 0, vec.Y, 0.0001)
	require.InDelta(t, -1, vec.Z, 0.0001)
}

func TestVec3_SetRotateZ(t *testing.T) {
	input := Vec3[float32]{X: 1, Y: 0, Z: 0}

	var vec Vec3[float32]
	vec.SetRotateZ(&input, math.Pi/2.0)

	require.InDelta(t, 0, vec.X, 0.0001)
	require.InDelta(t, 1, vec.Y, 0.0001)
	require.InDelta(t, 0, vec.Z, 0.0001)
}

func TestVec3_SetNormalizeVec3(t *testing.T) {
	input := Vec3[float32]{X: 5, Y: 0, Z: 0}

	var vec Vec3[float32]
	vec.SetNormalizeVec3(&input)

	require.InDelta(t, 1.0, vec.X, 0.0001)
	require.InDelta(t, 0, vec.Y, 0.0001)
	require.InDelta(t, 0, vec.Z, 0.0001)
}

func TestVec3_AddVec3(t *testing.T) {
	v1 := Vec3[float32]{X: 5, Y: 9, Z: 13}
	v2 := Vec3[float32]{X: 1, Y: 2, Z: 4}

	v1.AddVec3(&v2)

	require.InDelta(t, 6.0, v1.X, 0.0001)
	require.InDelta(t, 11.0, v1.Y, 0.0001)
	require.InDelta(t, 17.0, v1.Z, 0.0001)
}

func TestVec3_SetScale(t *testing.T) {
	input := Vec3[float32]{X: 1, Y: 2, Z: 4}

	var vec Vec3[float32]
	vec.SetScale(&input, 3.0)

	require.InDelta(t, 3.0, vec.X, 0.0001)
	require.InDelta(t, 6.0, vec.Y, 0.0001)
	require.InDelta(t, 12.0, vec.Z, 0.0001)
}

func TestVec3_Scale(t *testing.T) {
	vec := Vec3[float32]{X: 1, Y: 2, Z: 4}
	vec.Scale(3.0)

	require.InDelta(t, 3.0, vec.X, 0.0001)
	require.InDelta(t, 6.0, vec.Y, 0.0001)
	require.InDelta(t, 12.0, vec.Z, 0.0001)
}

func TestVec3_SetLerp(t *testing.T) {
	v1 := Vec3[float32]{X: 1, Y: 2, Z: 4}
	v2 := Vec3[float32]{X: 3, Y: 6, Z: 12}

	var vec Vec3[float32]
	vec.SetLerp(&v1, &v2, 0.5)

	require.InDelta(t, 2.0, vec.X, 0.0001)
	require.InDelta(t, 4.0, vec.Y, 0.0001)
	require.InDelta(t, 8.0, vec.Z, 0.0001)
}

func TestVec3_Lerp(t *testing.T) {
	v1 := Vec3[float32]{X: 1, Y: 2, Z: 4}
	v2 := Vec3[float32]{X: 3, Y: 6, Z: 12}

	v1.Lerp(&v2, 0.5)

	require.InDelta(t, 2.0, v1.X, 0.0001)
	require.InDelta(t, 4.0, v1.Y, 0.0001)
	require.InDelta(t, 8.0, v1.Z, 0.0001)
}

func TestVec3_SetClosestPointOnLine(t *testing.T) {
	point := Vec3[float32]{X: 1, Y: 0.9, Z: 0}
	endpoint1 := Vec3[float32]{X: 0, Y: 2.0, Z: 0}
	endpoint2 := Vec3[float32]{X: 0, Y: -2.0, Z: 0}

	var vec Vec3[float32]
	vec.SetClosestPointOnLine(&point, &endpoint1, &endpoint2)

	require.InDelta(t, 0, vec.X, 0.0001)
	require.InDelta(t, 0.9, vec.Y, 0.0001)
	require.InDelta(t, 0, vec.Z, 0.0001)
}

func TestVec3_Orthonormalize(t *testing.T) {
	v1 := Vec3[float32]{X: 3, Y: 0.15, Z: 0}
	v2 := Vec3[float32]{X: 0, Y: 1, Z: 0}

	v1.Orthonormalize(&v2)

	require.InDelta(t, 1.0, v1.Len(), 0.0001)
	require.InDelta(t, 0.0, v1.DotProduct(&v2), 0.0001)
}

func TestVec3_Rotate(t *testing.T) {
	v1 := Vec3[float32]{X: 1, Y: 0, Z: 0}
	axis := Vec3[float32]{X: 0, Y: 1, Z: 0}

	v1.Rotate(math.Pi/2.0, &axis)

	require.InDelta(t, 0, v1.X, 0.0001)
	require.InDelta(t, 0, v1.Y, 0.0001)
	require.InDelta(t, -1.0, v1.Z, 0.0001)
}

func TestVec3_RotateX(t *testing.T) {
	v1 := Vec3[float32]{X: 0, Y: 1, Z: 0}
	v1.RotateX(math.Pi / 2.0)

	require.InDelta(t, 0, v1.X, 0.0001)
	require.InDelta(t, 0, v1.Y, 0.0001)
	require.InDelta(t, 1, v1.Z, 0.0001)
}

func TestVec3_RotateY(t *testing.T) {
	v1 := Vec3[float32]{X: -1, Y: 0, Z: 0}
	v1.RotateY(math.Pi / 2.0)

	require.InDelta(t, 0, v1.X, 0.0001)
	require.InDelta(t, 0, v1.Y, 0.0001)
	require.InDelta(t, 1, v1.Z, 0.0001)
}

func TestVec3_RotateZ(t *testing.T) {
	v1 := Vec3[float32]{X: -1, Y: 0, Z: 0}
	v1.RotateZ(math.Pi / 2.0)

	require.InDelta(t, 0, v1.X, 0.0001)
	require.InDelta(t, -1, v1.Y, 0.0001)
	require.InDelta(t, 0, v1.Z, 0.0001)
}

func TestVec3_SetTriangleNormal(t *testing.T) {
	p1 := Vec3[float32]{X: 0, Y: 0, Z: 1}
	p2 := Vec3[float32]{X: 0, Y: 1, Z: 0}
	p3 := Vec3[float32]{X: 0, Y: 0, Z: -1}

	var vec Vec3[float32]
	vec.SetTriangleNormal(&p1, &p2, &p3)

	require.InDelta(t, -1, vec.X, 0.0001)
	require.InDelta(t, 0, vec.Y, 0.0001)
	require.InDelta(t, 0, vec.Z, 0.0001)
}
