package math

import "simd"

func TransformVec(vec *Vec4, mat *Mat4x4, out *Vec4) {
	mat0 := simd.LoadFloat32s(mat[0][:])
	x := simd.BroadcastFloat32s(vec.X())
	outVec := mat0.Mul(x)

	mat1 := simd.LoadFloat32s(mat[1][:])
	y := simd.BroadcastFloat32s(vec.Y())
	outVec = outVec.Add(mat1.Mul(y))

	mat2 := simd.LoadFloat32s(mat[2][:])
	z := simd.BroadcastFloat32s(vec.Z())
	outVec = outVec.Add(mat2.Mul(z))

	mat3 := simd.LoadFloat32s(mat[3][:])
	w := simd.BroadcastFloat32s(vec.W())
	outVec = outVec.Add(mat3.Mul(w))

	outVec.Store(out[:])
}

func MulMatrix(left *Mat4x4, right *Mat4x4, out *Mat4x4) {
	left0 := simd.LoadFloat32s(left[0][:])
	right0 := simd.BroadcastFloat32s(right[0][0])
	out0 := left0.Mul(right0)

	left1 := simd.LoadFloat32s(left[1][:])
	right1 := simd.BroadcastFloat32s(right[0][1])
	out0 = out0.Add(left1.Mul(right1))

	left2 := simd.LoadFloat32s(left[2][:])
	right2 := simd.BroadcastFloat32s(right[0][2])
	out0 = out0.Add(left2.Mul(right2))

	left3 := simd.LoadFloat32s(left[3][:])
	right3 := simd.BroadcastFloat32s(right[0][3])
	out0 = out0.Add(left3.Mul(right3))

	left0 = simd.LoadFloat32s(left[0][:])
	right0 = simd.BroadcastFloat32s(right[1][0])
	out1 := left0.Mul(right0)

	left1 = simd.LoadFloat32s(left[1][:])
	right1 = simd.BroadcastFloat32s(right[1][1])
	out1 = out1.Add(left1.Mul(right1))

	left2 = simd.LoadFloat32s(left[2][:])
	right2 = simd.BroadcastFloat32s(right[1][2])
	out1 = out1.Add(left2.Mul(right2))

	left3 = simd.LoadFloat32s(left[3][:])
	right3 = simd.BroadcastFloat32s(right[1][3])
	out1 = out1.Add(left3.Mul(right3))

	left0 = simd.LoadFloat32s(left[0][:])
	right0 = simd.BroadcastFloat32s(right[2][0])
	out2 := left0.Mul(right0)

	left1 = simd.LoadFloat32s(left[1][:])
	right1 = simd.BroadcastFloat32s(right[2][1])
	out2 = out2.Add(left1.Mul(right1))

	left2 = simd.LoadFloat32s(left[2][:])
	right2 = simd.BroadcastFloat32s(right[2][2])
	out2 = out2.Add(left2.Mul(right2))

	left3 = simd.LoadFloat32s(left[3][:])
	right3 = simd.BroadcastFloat32s(right[2][3])
	out2 = out2.Add(left3.Mul(right3))

	left0 = simd.LoadFloat32s(left[0][:])
	right0 = simd.BroadcastFloat32s(right[3][0])
	out3 := left0.Mul(right0)

	left1 = simd.LoadFloat32s(left[1][:])
	right1 = simd.BroadcastFloat32s(right[3][1])
	out3 = out3.Add(left1.Mul(right1))

	left2 = simd.LoadFloat32s(left[2][:])
	right2 = simd.BroadcastFloat32s(right[3][2])
	out3 = out3.Add(left2.Mul(right2))

	left3 = simd.LoadFloat32s(left[3][:])
	right3 = simd.BroadcastFloat32s(right[3][3])
	out3 = out3.Add(left3.Mul(right3))

	out0.Store(out[0][:])
	out1.Store(out[1][:])
	out2.Store(out[2][:])
	out3.Store(out[3][:])
}
