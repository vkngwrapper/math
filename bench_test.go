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
var g3nGlobalOutQuat math32.Quaternion
var g3nGlobalOutScalar float32

var mathglGlobalOutVec4 mgl32.Vec4
var mathglGlobalOutVec3 mgl32.Vec3
var mathglGlobalOutQuat mgl32.Quat
var mathglGlobalOutScalar float32

var go3dGlobalOutVec4 vec4.T
var go3dGlobalOutVec3 vec3.T
var go3dGlobalOutQuat quaternion.T
var go3dGlobalOutScalar float32

var vkngMathOutVec4 Vec4[float32]
var vkngMathOutVec3 Vec3[float32]
var vkngMathOutQuat Quaternion[float32]
var vkngMathOutScalar float32

func BenchmarkTransformVec4G3n(b *testing.B) {
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		transform := math32.Vector4{5, 10, 15, 1}

		var mat math32.Matrix4
		mat.Identity()
		mat.Multiply(mat.MakeTranslation(1, 1, 1)).
			Multiply(mat.MakeRotationY(1)).
			Multiply(mat.MakeScale(1.5, 1.5, 1.5))
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

func BenchmarkVkngMathTransform(b *testing.B) {
	b.ReportAllocs()

	for i := 0; i < b.N; i++ {
		translate := Vec3[float32]{1, 1, 1}
		scale := Vec3[float32]{1.5, 1.5, 1.5}

		var mat Mat4x4[float32]
		mat.SetIdentity().Translate(&translate).RotateY(1).Scale(&scale)

		transform := Vec4[float32]{5, 10, 15, 1}
		transform.Transform(&mat)

		vkngMathOutVec4 = transform
	}
}

func BenchmarkRotateVec3G3n(b *testing.B) {
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		transform := math32.Vector3{5, 10, 15}
		axis := math32.Vector3{1, 0, 0}

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
		transform := math32.Vector3{5, 10, 15}
		euler := math32.Vector3{math.Pi, 0, math.Pi / 2}

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
