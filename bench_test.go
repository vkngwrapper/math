package math

import (
	"github.com/ungerik/go3d/mat4"
	"github.com/ungerik/go3d/vec3"
	"testing"
)
import "github.com/g3n/engine/math32"
import "github.com/go-gl/mathgl/mgl32"

var g3nGlobalOut math32.Vector3

func BenchmarkG3nTransform(b *testing.B) {
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		transform := math32.Vector3{5, 10, 15}

		var mat math32.Matrix4
		mat.Identity()
		mat.Multiply(mat.MakeTranslation(1, 1, 1)).
			Multiply(mat.MakeRotationY(1)).
			Multiply(mat.MakeScale(1.5, 1.5, 1.5))
		g3nGlobalOut = *(transform.ApplyMatrix4(&mat))
	}
}

var mathglGlobalOut mgl32.Vec4

func BenchmarkMathGLTransform(b *testing.B) {
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {

		mat := mgl32.Ident4().
			Mul4(mgl32.Translate3D(1, 1, 1)).
			Mul4(mgl32.HomogRotate3DY(1)).
			Mul4(mgl32.Scale3D(1.5, 1.5, 1.5))
		transform := mgl32.Vec4{5.0, 10.0, 15.0, 1.0}
		transform = mat.Mul4x1(transform)

		mathglGlobalOut = transform
	}
}

var go3dGlobalOut vec3.T

func BenchmarkGo3DTransform(b *testing.B) {
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
		transform := vec3.T{5, 10, 15}
		mat.TransformVec3(&transform)

		go3dGlobalOut = transform
	}
}

var vkngMathOut Vec4[float32]

func BenchmarkVkngMathTransform(b *testing.B) {
	b.ReportAllocs()

	for i := 0; i < b.N; i++ {
		translate := Vec3[float32]{1, 1, 1}
		scale := Vec3[float32]{1.5, 1.5, 1.5}

		var mat Mat4x4[float32]
		mat.SetIdentity().Translate(&translate).RotateY(1).Scale(&scale)

		transform := Vec4[float32]{5, 10, 15, 1}
		transform.Transform(&mat)

		vkngMathOut = transform
	}
}
