package math

import (
	"github.com/stretchr/testify/require"
	"math"
	"testing"
)

func TestEqual4x4_True(t *testing.T) {
	mat1 := Mat4x4[float32]{
		{0.1, 0.2, 0.3, 0.4},
		{0.5, 0.6, 0.7, 0.8},
		{0.9, 1.0, 1.1, 1.2},
		{1.3, 1.4, 1.5, 1.6},
	}
	mat2 := Mat4x4[float32]{
		{0.1, 0.2, 0.3, 0.4},
		{0.5, 0.6, 0.7, 0.8},
		{0.9, 1.0, 1.1, 1.2},
		{1.3, 1.4, 1.5, 1.6},
	}
	require.True(t, mat1.Equal(&mat2, 0.0001))
}

func TestEqual4x4_False(t *testing.T) {
	mat1 := Mat4x4[float32]{
		{0.1, 0.2, 0.3, 0.4},
		{0.5, 0.6, 0.7, 0.8},
		{0.9, 1.0, 1.1, 1.2},
		{1.3, 1.4, 1.5, 1.6},
	}
	mat2 := Mat4x4[float32]{
		{0.1, 0.2, 0.3, 0.4},
		{0.5, 0.6, 0.7, 0.8},
		{0.9, 1.0, 1.2, 1.2},
		{1.3, 1.4, 1.5, 1.6},
	}
	require.False(t, mat1.Equal(&mat2, 0.0001))
}

func TestInverse4x4(t *testing.T) {
	mat1 := Mat4x4[float32]{
		{0.6, 0.2, 0.3, 0.4},
		{0.2, 0.7, 0.5, 0.3},
		{0.3, 0.5, 0.7, 0.2},
		{0.4, 0.3, 0.2, 0.6},
	}
	var inverse Mat4x4[float32]
	inverse.SetMat4x4(&mat1).Inverse()
	var identity Mat4x4[float32]
	identity.SetIdentity()

	require.True(t, mat1.MultMat4x4(&inverse).Equal(&identity, 0.0001))
}

func TestTranspose4x4(t *testing.T) {
	mat := Mat4x4[float32]{
		{0, 1, 2, 3},
		{4, 5, 6, 7},
		{8, 9, 10, 11},
		{12, 13, 14, 15},
	}
	expected := Mat4x4[float32]{
		{0, 4, 8, 12},
		{1, 5, 9, 13},
		{2, 6, 10, 14},
		{3, 7, 11, 15},
	}

	require.True(t, expected.Equal(mat.Transpose(), 0.0001))
}

func TestMat4x4_SetFrustum(t *testing.T) {
	var frustum Mat4x4[float32]
	frustum.SetFrustum(3, 8, 6, 9, 5, 7)

	point := Vec4[float32]{5.5, 7.5, -5, 1}
	point.Transform(&frustum)

	require.InDelta(t, 0, point.X/point.W, 0.0001)
	require.InDelta(t, 0, point.Y/point.W, 0.0001)
	require.InDelta(t, 0, point.Z/point.W, 0.0001)

	frustum.SetFrustum(-6, -3, -5, -3, 3, 5)
	point = Vec4[float32]{-3, -3, -3, 1}
	point.Transform(&frustum)

	require.InDelta(t, 1, point.X/point.W, 0.0001)
	require.InDelta(t, 1, point.Y/point.W, 0.0001)
	require.InDelta(t, 0, point.Z/point.W, 0.0001)
}

func TestMat4x4_SetPerspective(t *testing.T) {
	var perspective Mat4x4[float32]
	perspective.SetPerspective(math.Pi*0.25, 2, 3, 5)

	point := Vec4[float32]{-2.4853, 1.2426, -3, 1}
	point.Transform(&perspective)

	require.InDelta(t, -1, point.X/point.W, 0.0001)
	require.InDelta(t, -1, point.Y/point.W, 0.0001)
	require.InDelta(t, 0, point.Z/point.W, 0.0001)

	point = Vec4[float32]{4.142, -2.071, -5, 1}
	point.Transform(&perspective)

	require.InDelta(t, 1, point.X/point.W, 0.0001)
	require.InDelta(t, 1, point.Y/point.W, 0.0001)
	require.InDelta(t, 1, point.Z/point.W, 0.0001)
}

func TestMat4x4_SetLookAt(t *testing.T) {
	var viewMat Mat4x4[float32]
	viewMat.SetLookAt(&Vec3[float32]{1, 5, 3}, &Vec3[float32]{1, 0, 3}, &Vec3[float32]{1, 0, 0})

	point := Vec4[float32]{1, 0, 3, 1}
	point.Transform(&viewMat)

	require.InDelta(t, float32(0), point.X, 0.0001)
	require.InDelta(t, float32(0), point.Y, 0.0001)
	require.InDelta(t, float32(-5), point.Z, 0.0001)
	require.InDelta(t, float32(1), point.W, 0.0001)
}

func TestMat4x4_SetInfinitePerspective(t *testing.T) {
	var proj Mat4x4[float32]
	proj.SetInfinitePerspective(math.Pi*0.25, 2, 3)

	point := Vec4[float32]{0, 0, -500, 0}
	point.Transform(&proj)

	require.InDelta(t, 0, point.X/point.W, 0.0001)
	require.InDelta(t, 0, point.Y/point.W, 0.0001)
	require.InDelta(t, 0.9999, point.Z/point.W, 0.0001)
}

func TestMat4x4_SetOrthographic(t *testing.T) {
	var proj Mat4x4[float32]
	proj.SetOrthographic(-100, 500, -50, 350, 3, 100)

	point := Vec4[float32]{200, 150, -100, 1}
	point.Transform(&proj)

	require.InDelta(t, 0, point.X/point.W, 0.0001)
	require.InDelta(t, 0, point.Y/point.W, 0.0001)
	require.InDelta(t, 1, point.Z/point.W, 0.0001)
}

func TestMat4x4_SetOrthographic2D(t *testing.T) {
	var proj Mat4x4[float32]
	proj.SetOrthographic2D(-100, 500, -50, 350)

	point := Vec4[float32]{-100, -50, 100, 1}
	point.Transform(&proj)

	require.InDelta(t, -1, point.X/point.W, 0.0001)
	require.InDelta(t, -1, point.Y/point.W, 0.0001)
	require.InDelta(t, 100, point.Z/point.W, 0.0001)
}

func TestMat4x4_SetAxisAngle(t *testing.T) {
	axis := Vec3[float32]{0, 1, 0}
	var rotation Mat4x4[float32]
	rotation.SetRotateAroundAxis(&axis, math.Pi)

	point := Vec4[float32]{1, 0, 1, 1}
	point.Transform(&rotation)

	require.InDelta(t, -1, point.X/point.W, 0.0001)
	require.InDelta(t, 0, point.Y/point.W, 0.0001)
	require.InDelta(t, -1, point.Z/point.W, 0.0001)
}

func TestMat4x4_SetAffineMat4x4(t *testing.T) {
	axis := Vec3[float32]{0, 0, 1}
	var rotation Mat4x4[float32]
	rotation.SetAxisAngle(&axis, math.Pi/2.0)

	var rotation2 Mat4x4[float32]
	rotation2.SetAffineMat4x4(&rotation)

	var outAxis Vec3[float32]
	var angle float64
	rotation2.GetAxisAngle(&outAxis, &angle)
	require.Equal(t, float32(0.0), axis.X)
	require.Equal(t, float32(0.0), axis.Y)
	require.Equal(t, float32(1.0), axis.Z)

	require.InDelta(t, math.Pi/2.0, angle, 0.0001)
}

func TestMat4x4_SetOrientation(t *testing.T) {
	target := Vec3[float32]{1, 0, 0}
	source := Vec3[float32]{0, 0, 1}

	var rotation Mat4x4[float32]
	rotation.SetOrientation(&source, &target)

	point := Vec4[float32]{1, 0, 1, 1}
	point.Transform(&rotation)

	require.InDelta(t, 1, point.X/point.W, 0.0001)
	require.InDelta(t, 0, point.Y/point.W, 0.0001)
	require.InDelta(t, -1, point.Z/point.W, 0.0001)
}

func TestMat4x4_SetTranslation(t *testing.T) {
	var translation Mat4x4[float32]
	translation.SetTranslation(1, 2, 3)

	point := Vec4[float32]{5, 6, 7, 1}
	point.Transform(&translation)

	require.InDelta(t, 6, point.X/point.W, 0.0001)
	require.InDelta(t, 8, point.Y/point.W, 0.0001)
	require.InDelta(t, 10, point.Z/point.W, 0.0001)
}

func TestMat4x4_Translate_Identity(t *testing.T) {
	var translation Mat4x4[float32]
	translation.SetIdentity()
	translation.Translate(1, 2, 3)

	point := Vec4[float32]{5, 6, 7, 1}
	point.Transform(&translation)

	require.InDelta(t, 6, point.X/point.W, 0.0001)
	require.InDelta(t, 8, point.Y/point.W, 0.0001)
	require.InDelta(t, 10, point.Z/point.W, 0.0001)
}

func TestMat4x4_Translate_FromScale(t *testing.T) {
	var translation Mat4x4[float32]
	translation.SetScale(2, 3, 2)
	translation.Translate(1, 2, 3)

	point := Vec4[float32]{5, 6, 7, 1}
	point.Transform(&translation)

	require.InDelta(t, 11, point.X/point.W, 0.0001)
	require.InDelta(t, 20, point.Y/point.W, 0.0001)
	require.InDelta(t, 17, point.Z/point.W, 0.0001)
}

func TestMat4x4_Scale_Identity(t *testing.T) {
	var scale Mat4x4[float32]
	scale.SetIdentity()
	scale.Scale(2, 3, 4)

	point := Vec4[float32]{3, 5, 7, 1}
	point.Transform(&scale)

	require.InDelta(t, 6, point.X/point.W, 0.0001)
	require.InDelta(t, 15, point.Y/point.W, 0.0001)
	require.InDelta(t, 28, point.Z/point.W, 0.0001)
}

func TestMat4x4_Scale_FromTranslate(t *testing.T) {
	var scale Mat4x4[float32]
	scale.SetTranslation(1, 2, 3)
	scale.Scale(2, 3, 4)

	point := Vec4[float32]{3, 5, 7, 1}
	point.Transform(&scale)

	require.InDelta(t, 8, point.X/point.W, 0.0001)
	require.InDelta(t, 21, point.Y/point.W, 0.0001)
	require.InDelta(t, 40, point.Z/point.W, 0.0001)
}

func TestMat4x4_SetScale(t *testing.T) {
	var scale Mat4x4[float32]
	scale.SetScale(1.5, 2.5, 3.5)

	point := Vec4[float32]{1, 2, 3, 1}
	point.Transform(&scale)

	require.InDelta(t, 1.5, point.X/point.W, 0.0001)
	require.InDelta(t, 5, point.Y/point.W, 0.0001)
	require.InDelta(t, 10.5, point.Z/point.W, 0.0001)
}

func TestMat4x4_SetShearX(t *testing.T) {
	var shear Mat4x4[float32]
	shear.SetShearX(0.5, 1)

	point := Vec4[float32]{1, 0, 0, 1}
	point.Transform(&shear)

	require.InDelta(t, 1.0, point.X/point.W, 0.0001)
	require.InDelta(t, 0.5, point.Y/point.W, 0.0001)
	require.InDelta(t, 1.0, point.Z/point.W, 0.0001)
}

func TestMat4x4_ShearX_Identity(t *testing.T) {
	var shear Mat4x4[float32]
	shear.SetIdentity()
	shear.ShearX(0.5, 1)

	point := Vec4[float32]{1, 0, 0, 1}
	point.Transform(&shear)

	require.InDelta(t, 1.0, point.X/point.W, 0.0001)
	require.InDelta(t, 0.5, point.Y/point.W, 0.0001)
	require.InDelta(t, 1.0, point.Z/point.W, 0.0001)
}

func TestMat4x4_ShearX_FromTranslate(t *testing.T) {
	var shear Mat4x4[float32]
	shear.SetTranslation(3, 0, 0)
	shear.ShearX(0.5, 1)

	point := Vec4[float32]{1, 0, 0, 1}
	point.Transform(&shear)

	require.InDelta(t, 4.0, point.X/point.W, 0.0001)
	require.InDelta(t, 2.0, point.Y/point.W, 0.0001)
	require.InDelta(t, 4.0, point.Z/point.W, 0.0001)
}

func TestMat4x4_SetShearY(t *testing.T) {
	var shear Mat4x4[float32]
	shear.SetShearY(0.5, 1)

	point := Vec4[float32]{0, 1, 0, 1}
	point.Transform(&shear)

	require.InDelta(t, 0.5, point.X/point.W, 0.0001)
	require.InDelta(t, 1.0, point.Y/point.W, 0.0001)
	require.InDelta(t, 1.0, point.Z/point.W, 0.0001)
}

func TestMat4x4_ShearY_Identity(t *testing.T) {
	var shear Mat4x4[float32]
	shear.SetIdentity()
	shear.ShearY(0.5, 1)

	point := Vec4[float32]{0, 1, 0, 1}
	point.Transform(&shear)

	require.InDelta(t, 0.5, point.X/point.W, 0.0001)
	require.InDelta(t, 1.0, point.Y/point.W, 0.0001)
	require.InDelta(t, 1.0, point.Z/point.W, 0.0001)
}

func TestMat4x4_ShearY_FromTranslate(t *testing.T) {
	var shear Mat4x4[float32]
	shear.SetTranslation(0, 3, 0)
	shear.ShearY(0.5, 1)

	point := Vec4[float32]{0, 1, 0, 1}
	point.Transform(&shear)

	require.InDelta(t, 2.0, point.X/point.W, 0.0001)
	require.InDelta(t, 4.0, point.Y/point.W, 0.0001)
	require.InDelta(t, 4.0, point.Z/point.W, 0.0001)
}

func TestMat4x4_SetShearZ(t *testing.T) {
	var shear Mat4x4[float32]
	shear.SetShearZ(0.5, 1)

	point := Vec4[float32]{0, 0, 1, 1}
	point.Transform(&shear)

	require.InDelta(t, 0.5, point.X/point.W, 0.0001)
	require.InDelta(t, 1.0, point.Y/point.W, 0.0001)
	require.InDelta(t, 1.0, point.Z/point.W, 0.0001)
}

func TestMat4x4_ShearZ_Identity(t *testing.T) {
	var shear Mat4x4[float32]
	shear.SetIdentity()
	shear.ShearZ(0.5, 1)

	point := Vec4[float32]{0, 0, 1, 1}
	point.Transform(&shear)

	require.InDelta(t, 0.5, point.X/point.W, 0.0001)
	require.InDelta(t, 1.0, point.Y/point.W, 0.0001)
	require.InDelta(t, 1.0, point.Z/point.W, 0.0001)
}

func TestMat4x4_ShearZ_FromTranslate(t *testing.T) {
	var shear Mat4x4[float32]
	shear.SetTranslation(1, 1, 3)
	shear.ShearZ(0.5, 1)

	point := Vec4[float32]{0, 0, 1, 1}
	point.Transform(&shear)

	require.InDelta(t, 3.0, point.X/point.W, 0.0001)
	require.InDelta(t, 5.0, point.Y/point.W, 0.0001)
	require.InDelta(t, 4.0, point.Z/point.W, 0.0001)
}

func TestMat4x4_SetRotationY(t *testing.T) {
	var rotation Mat4x4[float32]
	rotation.SetRotationY(math.Pi / 2)

	point := Vec4[float32]{1, 0, 0, 1}
	point.Transform(&rotation)

	require.InDelta(t, 0, point.X/point.W, 0.0001)
	require.InDelta(t, 0, point.Y/point.W, 0.0001)
	require.InDelta(t, -1, point.Z/point.W, 0.0001)
}

func TestMat4x4_RotateY_Identity(t *testing.T) {
	var rotation Mat4x4[float32]
	rotation.SetIdentity()
	rotation.RotateY(math.Pi / 2)

	point := Vec4[float32]{1, 0, 0, 1}
	point.Transform(&rotation)

	require.InDelta(t, 0, point.X/point.W, 0.0001)
	require.InDelta(t, 0, point.Y/point.W, 0.0001)
	require.InDelta(t, -1, point.Z/point.W, 0.0001)
}

func TestMat4x4_RotateY_FromTranslate(t *testing.T) {
	var rotation Mat4x4[float32]
	rotation.SetIdentity()
	rotation.SetTranslation(1, 0, 0)
	rotation.RotateY(math.Pi / 2)

	point := Vec4[float32]{0, 0, 0, 1}
	point.Transform(&rotation)

	require.InDelta(t, 0, point.X/point.W, 0.0001)
	require.InDelta(t, 0, point.Y/point.W, 0.0001)
	require.InDelta(t, -1, point.Z/point.W, 0.0001)
}

func TestMat4x4_SetRotationX(t *testing.T) {
	var rotation Mat4x4[float32]
	rotation.SetRotationX(math.Pi / 2)

	point := Vec4[float32]{0, 1, 0, 1}
	point.Transform(&rotation)

	require.InDelta(t, 0, point.X/point.W, 0.0001)
	require.InDelta(t, 0, point.Y/point.W, 0.0001)
	require.InDelta(t, 1, point.Z/point.W, 0.0001)
}

func TestMat4x4_RotateX_Identity(t *testing.T) {
	var rotation Mat4x4[float32]
	rotation.SetIdentity()
	rotation.RotateX(math.Pi / 2)

	point := Vec4[float32]{0, 1, 0, 1}
	point.Transform(&rotation)

	require.InDelta(t, 0, point.X/point.W, 0.0001)
	require.InDelta(t, 0, point.Y/point.W, 0.0001)
	require.InDelta(t, 1, point.Z/point.W, 0.0001)
}

func TestMat4x4_RotateX_FromTranslate(t *testing.T) {
	var rotation Mat4x4[float32]
	rotation.SetTranslation(0, 1, 0)
	rotation.RotateX(math.Pi / 2)

	point := Vec4[float32]{0, 0, 0, 1}
	point.Transform(&rotation)

	require.InDelta(t, 0, point.X/point.W, 0.0001)
	require.InDelta(t, 0, point.Y/point.W, 0.0001)
	require.InDelta(t, 1, point.Z/point.W, 0.0001)
}

func TestMat4x4_SetRotationZ(t *testing.T) {
	var rotation Mat4x4[float32]
	rotation.SetRotationZ(math.Pi / 2)

	point := Vec4[float32]{1, 1, 0, 1}
	point.Transform(&rotation)

	require.InDelta(t, -1, point.X/point.W, 0.0001)
	require.InDelta(t, 1, point.Y/point.W, 0.0001)
	require.InDelta(t, 0, point.Z/point.W, 0.0001)
}

func TestMat4x4_RotateZ_Identity(t *testing.T) {
	var rotation Mat4x4[float32]
	rotation.SetIdentity()
	rotation.RotateZ(math.Pi / 2)

	point := Vec4[float32]{1, 1, 0, 1}
	point.Transform(&rotation)

	require.InDelta(t, -1, point.X/point.W, 0.0001)
	require.InDelta(t, 1, point.Y/point.W, 0.0001)
	require.InDelta(t, 0, point.Z/point.W, 0.0001)
}

func TestMat4x4_RotateZ_FromTranslate(t *testing.T) {
	var rotation Mat4x4[float32]
	rotation.SetTranslation(1, 1, 0)
	rotation.RotateZ(math.Pi / 2)

	point := Vec4[float32]{0, 0, 0, 1}
	point.Transform(&rotation)

	require.InDelta(t, -1, point.X/point.W, 0.0001)
	require.InDelta(t, 1, point.Y/point.W, 0.0001)
	require.InDelta(t, 0, point.Z/point.W, 0.0001)
}

func TestMat4x4_SetRotationEulers(t *testing.T) {
	var eulers Mat4x4[float32]
	eulers.SetRotationEulers(math.Pi/2, -math.Pi/2, math.Pi)

	point := Vec4[float32]{-1, 0, 0, 1}
	point.Transform(&eulers)

	require.InDelta(t, 0, point.X/point.W, 0.0001)
	require.InDelta(t, -1, point.Y/point.W, 0.0001)
	require.InDelta(t, 0, point.Z/point.W, 0.0001)
}

func TestMat4x4_GetAxisAngle(t *testing.T) {
	axis := Vec3[float32]{0, 1, 0}
	var axisAngleMat Mat4x4[float32]
	axisAngleMat.SetRotateAroundAxis(&axis, math.Pi/2)

	var eulerMat Mat4x4[float32]
	eulerMat.SetRotationEulers(math.Pi/2, 0, 0)
	require.True(t, eulerMat.Equal(&axisAngleMat, 0.0001))

	var angle float64

	axisAngleMat.GetAxisAngle(&axis, &angle)

	require.InDelta(t, 0, axis.X, 0.0001)
	require.InDelta(t, 1, axis.Y, 0.0001)
	require.InDelta(t, 0, axis.Z, 0.0001)
	require.InDelta(t, math.Pi/2, angle, 0.0001)
}

func TestMat4x4_SetRotationAroundAxis(t *testing.T) {
	axis := Vec3[float32]{0, 1, 0}
	var rotation Mat4x4[float32]
	rotation.SetRotateAroundAxis(&axis, math.Pi)

	point := Vec4[float32]{1, 0, 1, 1}
	point.Transform(&rotation)

	require.InDelta(t, -1, point.X/point.W, 0.0001)
	require.InDelta(t, 0, point.Y/point.W, 0.0001)
	require.InDelta(t, -1, point.Z/point.W, 0.0001)
}

func TestMat4x4_RotateAroundAxis_Identity(t *testing.T) {
	axis := Vec3[float32]{0, 1, 0}
	var rotation Mat4x4[float32]
	rotation.SetIdentity()
	rotation.RotateAroundAxis(&axis, math.Pi)

	point := Vec4[float32]{1, 0, 1, 1}
	point.Transform(&rotation)

	require.InDelta(t, -1, point.X/point.W, 0.0001)
	require.InDelta(t, 0, point.Y/point.W, 0.0001)
	require.InDelta(t, -1, point.Z/point.W, 0.0001)
}

func TestMat4x4_RotateAroundAxis_Translate(t *testing.T) {
	axis := Vec3[float32]{0, 1, 0}
	var rotation Mat4x4[float32]
	rotation.SetTranslation(1, 0, 1)
	rotation.RotateAroundAxis(&axis, math.Pi)

	point := Vec4[float32]{0, 0, 0, 1}
	point.Transform(&rotation)

	require.InDelta(t, -1, point.X/point.W, 0.0001)
	require.InDelta(t, 0, point.Y/point.W, 0.0001)
	require.InDelta(t, -1, point.Z/point.W, 0.0001)
}

func TestMat4x4_MultMatrix4x4(t *testing.T) {
	var mat1 Mat4x4[float32]
	mat1.SetRotationY(math.Pi / 2.0)

	var mat2 Mat4x4[float32]
	mat2.SetRotationY(math.Pi / 2.0)

	var expected Mat4x4[float32]
	expected.SetRotationY(math.Pi)

	mat1.MultMat4x4(&mat2)

	require.True(t, mat1.Equal(&expected, 0.0001))
	require.False(t, mat2.Equal(&expected, 0.0001))
}

func TestMat4x4_SetMultMatrix4x4(t *testing.T) {
	var mat1 Mat4x4[float32]
	mat1.SetRotationY(math.Pi / 2.0)

	var mat2 Mat4x4[float32]
	mat2.SetRotationY(math.Pi / 2.0)

	var expected Mat4x4[float32]
	expected.SetRotationY(math.Pi)

	var result Mat4x4[float32]
	result.SetMultMat4x4(&mat1, &mat2)

	require.True(t, result.Equal(&expected, 0.0001))
	require.False(t, mat1.Equal(&expected, 0.0001))
	require.False(t, mat2.Equal(&expected, 0.0001))
}

func TestMat4x4_InterpolateMatrix(t *testing.T) {
	var mat1 Mat4x4[float32]
	mat1.SetIdentity()

	var mat2 Mat4x4[float32]
	mat2.SetRotationY(math.Pi)

	var expected Mat4x4[float32]
	expected.SetRotationY(math.Pi / 2.0)

	mat1.InterpolateMat4x4(&mat2, 0.5)

	require.True(t, mat1.Equal(&expected, 0.0001))
	require.False(t, mat2.Equal(&expected, 0.0001))
}

func TestMat4x4_SetInterpolate(t *testing.T) {
	var mat1 Mat4x4[float32]
	mat1.SetIdentity()

	var mat2 Mat4x4[float32]
	mat2.SetRotationY(math.Pi)

	var expected Mat4x4[float32]
	expected.SetRotationY(math.Pi / 2.0)

	var result Mat4x4[float32]
	result.SetInterpolateMat4x4(&mat1, &mat2, 0.5)

	require.True(t, result.Equal(&expected, 0.0001))
	require.False(t, mat1.Equal(&expected, 0.0001))
	require.False(t, mat2.Equal(&expected, 0.0001))
}

func TestMat4x4_IsNormalized(t *testing.T) {
	var mat1 Mat4x4[float32]
	mat1.SetIdentity()

	require.True(t, mat1.IsNormalized(0.0001))

	mat1.Translate(5, 10, 15)
	require.False(t, mat1.IsNormalized(0.0001))
}

func TestMat4x4_ApplyTransform(t *testing.T) {
	var mat1 Mat4x4[float32]
	mat1.SetRotationY(math.Pi / 2.0)

	var mat2 Mat4x4[float32]
	mat2.SetRotationY(math.Pi / 2.0)

	var expected Mat4x4[float32]
	expected.SetRotationY(math.Pi)

	mat1.ApplyTransform(&mat2)

	require.True(t, mat1.Equal(&expected, 0.0001))
	require.False(t, mat2.Equal(&expected, 0.0001))
}
