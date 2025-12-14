package math

import (
	"math"
	"testing"
)

var vkngMathOutVec4 Vec4[float32]
var vkngMathOutVec3 Vec3[float32]
var vkngMathOutQuat Quaternion[float32]

func BenchmarkTransformVec4VkngMath(b *testing.B) {
	b.ReportAllocs()

	for i := 0; i < b.N; i++ {
		var mat Mat4x4[float32]
		mat.SetIdentity()
		mat.Translate(1, 1, 1)
		mat.RotateY(1)
		mat.Scale(1.5, 1.5, 1.5)

		transform := Vec4[float32]{5, 10, 15, 1}
		transform.Transform(&mat)

		vkngMathOutVec4 = transform
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

func BenchmarkRotateQuaternionVkngMath(b *testing.B) {
	b.ReportAllocs()

	for i := 0; i < b.N; i++ {
		transform := Vec3[float32]{5, 10, 15}
		var quat Quaternion[float32]
		quat.SetRotationEulers(0, math.Pi, math.Pi/2)
		transform.RotateWithQuaternion(&quat)

		vkngMathOutQuat = quat
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
		mat.SetMultMat4x4(&translateMat, &rotateMat)
		mat.MultMat4x4(&scaleMat)

		transform := Vec4[float32]{5, 10, 15, 1}
		transform.Transform(&mat)

		vkngMathOutVec4 = transform
	}
}
