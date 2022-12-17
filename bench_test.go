package math

import (
	"github.com/g3n/engine/math32"
	"github.com/go-gl/mathgl/mgl32"
	"github.com/ungerik/go3d/mat4"
	"github.com/ungerik/go3d/quaternion"
	"github.com/ungerik/go3d/vec3"
	"github.com/ungerik/go3d/vec4"
	"math"
	"testing"
)

var g3nGlobalOutVec4 math32.Vector4
var g3nGlobalOutVec3 math32.Vector3

var mathglGlobalOutVec4 mgl32.Vec4
var mathglGlobalOutVec3 mgl32.Vec3

var go3dGlobalOutVec4 vec4.T
var go3dGlobalOutVec3 vec3.T

var vkngMathOutVec4 Vec4[float32]
var vkngMathOutVec3 Vec3[float32]
var vkngMathOutQuat Quaternion[float32]

func BenchmarkTransformVec4G3n(b *testing.B) {
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		transform := math32.Vector4{X: 5, Y: 10, Z: 15, W: 1}

		var translateMat, rotateMat, scaleMat math32.Matrix4
		translateMat.MakeTranslation(1, 1, 1)
		rotateMat.MakeRotationY(1)
		scaleMat.MakeScale(1.5, 1.5, 1.5)

		var mat math32.Matrix4
		mat.MultiplyMatrices(
			&translateMat,
			&rotateMat,
		).Multiply(&scaleMat)

		g3nGlobalOutVec4 = *(transform.ApplyMatrix4(&mat))
	}
}

func BenchmarkTransformVec4MathGL(b *testing.B) {
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {

		mat := mgl32.Ident4().
			Mul4(mgl32.Translate3D(1, 1, 1)).
			Mul4(mgl32.HomogRotate3DY(1)).
			Mul4(mgl32.Scale3D(1.5, 1.5, 1.5))
		transform := mgl32.Vec4{5.0, 10.0, 15.0, 1.0}
		transform = mat.Mul4x1(transform)

		mathglGlobalOutVec4 = transform
	}
}

func BenchmarkTransformVec4Go3D(b *testing.B) {
	b.ReportAllocs()

	for i := 0; i < b.N; i++ {
		translate := vec3.T{1, 1, 1}
		scale := vec3.T{1.5, 1.5, 1.5}

		rotate := mat4.Ident
		rotate.AssignYRotation(1)

		mat := mat4.Ident
		mat.Translate(&translate).
			MultMatrix(&rotate).
			ScaleVec3(&scale)
		transform := vec4.T{5, 10, 15, 1}
		mat.TransformVec4(&transform)

		go3dGlobalOutVec4 = transform
	}
}

func BenchmarkTransformVec4VkngMath(b *testing.B) {
	b.ReportAllocs()

	for i := 0; i < b.N; i++ {
		var mat Mat4x4[float32]
		mat.SetIdentity().Translate(1, 1, 1).RotateY(1).Scale(1.5, 1.5, 1.5)

		transform := Vec4[float32]{5, 10, 15, 1}
		transform.Transform(&mat)

		vkngMathOutVec4 = transform
	}
}

func BenchmarkRotateVec3G3n(b *testing.B) {
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		transform := math32.Vector3{X: 5, Y: 10, Z: 15}
		axis := math32.Vector3{X: 1, Y: 0, Z: 0}

		transform.ApplyAxisAngle(&axis, math.Pi)

		g3nGlobalOutVec3 = transform
	}
}

func BenchmarkRotateVec3MathGL(b *testing.B) {
	b.ReportAllocs()

	for i := 0; i < b.N; i++ {
		transform := mgl32.Vec3{5, 10, 15}
		axis := mgl32.Vec3{1, 0, 0}
		quaternion := mgl32.QuatRotate(math.Pi, axis)

		mathglGlobalOutVec3 = quaternion.Rotate(transform)
	}
}

func BenchmarkRotateVec3Go3D(b *testing.B) {
	b.ReportAllocs()

	for i := 0; i < b.N; i++ {
		transform := vec3.T{5, 10, 15}
		quaternion := quaternion.FromXAxisAngle(math.Pi)

		quaternion.RotateVec3(&transform)
		go3dGlobalOutVec3 = transform
	}
}

func BenchmarkRotateVec3VkngMath(b *testing.B) {
	b.ReportAllocs()

	for i := 0; i < b.N; i++ {
		transform := Vec3[float32]{5, 10, 15}
		axis := Vec3[float32]{1, 0, 0}

		transform.Rotate(math.Pi, &axis)
		vkngMathOutVec3 = transform
	}
}

func BenchmarkRotateQuaternionG3n(b *testing.B) {
	b.ReportAllocs()

	for i := 0; i < b.N; i++ {
		transform := math32.Vector3{X: 5, Y: 10, Z: 15}
		euler := math32.Vector3{X: math.Pi, Y: 0, Z: math.Pi / 2}

		var quaternion math32.Quaternion

		quaternion.SetFromEuler(&euler)

		transform.ApplyQuaternion(&quaternion)

		g3nGlobalOutVec3 = transform
	}
}

func BenchmarkRotateQuaternionMathGL(b *testing.B) {
	b.ReportAllocs()

	for i := 0; i < b.N; i++ {
		transform := mgl32.Vec3{5, 10, 15}

		quaternion := mgl32.AnglesToQuat(math.Pi, 0, math.Pi/2, mgl32.XYZ)

		mathglGlobalOutVec3 = quaternion.Rotate(transform)
	}
}

func BenchmarkRotateQuaternionGo3D(b *testing.B) {
	b.ReportAllocs()

	for i := 0; i < b.N; i++ {
		transform := vec3.T{5, 10, 15}
		quaternion := quaternion.FromEulerAngles(0, math.Pi, math.Pi/2)

		quaternion.RotateVec3(&transform)
		go3dGlobalOutVec3 = transform
	}
}

func BenchmarkRotateQuaternionVkngMath(b *testing.B) {
	b.ReportAllocs()

	for i := 0; i < b.N; i++ {
		transform := Vec3[float32]{5, 10, 15}
		var quat Quaternion[float32]
		quat.SetEulerAngles(0, math.Pi, math.Pi/2)
		transform.RotateWithQuaternion(&quat)

		vkngMathOutQuat = quat
	}
}

func BenchmarkMatrixMultTransformG3n(b *testing.B) {
	b.ReportAllocs()

	for i := 0; i < b.N; i++ {
		transform := math32.Vector4{X: 5, Y: 10, Z: 15, W: 1}

		var translateMat, rotateMat, scaleMat math32.Matrix4
		translateMat.MakeTranslation(1, 1, 1)
		rotateMat.MakeRotationY(1)
		scaleMat.MakeScale(1.5, 1.5, 1.5)

		var mat math32.Matrix4
		mat.MultiplyMatrices(
			&translateMat,
			&rotateMat,
		).Multiply(&scaleMat)

		g3nGlobalOutVec4 = *(transform.ApplyMatrix4(&mat))
	}
}

func BenchmarkMatrixMultTransformMatGL(b *testing.B) {
	b.ReportAllocs()

	for i := 0; i < b.N; i++ {
		mat := mgl32.Ident4().
			Mul4(mgl32.Translate3D(1, 1, 1)).
			Mul4(mgl32.HomogRotate3DY(1)).
			Mul4(mgl32.Scale3D(1.5, 1.5, 1.5))
		transform := mgl32.Vec4{5.0, 10.0, 15.0, 1.0}
		transform = mat.Mul4x1(transform)

		mathglGlobalOutVec4 = transform
	}
}

func BenchmarkMatrixMultTransformGo3D(b *testing.B) {
	b.ReportAllocs()

	for i := 0; i < b.N; i++ {
		translate := vec3.T{1, 1, 1}
		scale := vec3.T{1.5, 1.5, 1.5}

		rotate := mat4.Ident
		rotate.AssignYRotation(1)

		translateMat := mat4.Ident
		translateMat.Translate(&translate)

		scaleMat := mat4.Ident
		scaleMat.ScaleVec3(&scale)

		var mat mat4.T
		mat.AssignMul(&translateMat, &rotate).
			MultMatrix(&scaleMat)
		transform := vec4.T{5, 10, 15, 1}
		mat.TransformVec4(&transform)

		go3dGlobalOutVec4 = transform
	}
}

func BenchmarkMatrixMultTransformVkngMath(b *testing.B) {
	b.ReportAllocs()

	for i := 0; i < b.N; i++ {

		var translateMat, rotateMat, scaleMat Mat4x4[float32]
		translateMat.SetTranslation(1, 1, 1)
		rotateMat.SetRotationY(1)
		scaleMat.SetScale(1.5, 1.5, 1.5)

		var mat Mat4x4[float32]
		mat.SetMultMatrix4x4(&translateMat, &rotateMat).MultMatrix4x4(&scaleMat)

		transform := Vec4[float32]{5, 10, 15, 1}
		transform.Transform(&mat)

		vkngMathOutVec4 = transform
	}
}
