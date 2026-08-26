package math

import "math"

// Vec4 is a three-element vector of floating point compatible values that can be used for 3d
// vector arithmetic and transforms
type Vec4 [4]float32

func (v *Vec4) X() float32 {
	return v[0]
}

func (v *Vec4) Y() float32 {
	return v[1]
}

func (v *Vec4) Z() float32 {
	return v[2]
}

func (v *Vec4) W() float32 {
	return v[3]
}

// SetVec4 overwrites the contents of this vector with the contents of a 4-element vector
//
// in - The vector to initialize from
func (v *Vec4) SetVec4(in *Vec4) {
	v[0] = in[0]
	v[1] = in[1]
	v[2] = in[2]
	v[3] = in[3]
}

// SetVec3 overwrites the first three elements of this vector with the first three elements of a 3-element vector.
// The fourth element of this vector is set to 1.
//
// in - The vector to initialize from
func (v *Vec4) SetVec3(in *Vec3) {
	v[0] = in.X
	v[1] = in.Y
	v[2] = in.Z
	v[3] = 1
}

// SetVec2 overwrites the first two elements of this vector with the first two elements of a 2-element vector.
// The third element of this vector is set to 0 and the fourth is set to 1.
//
// in - The vector to initialize from
func (v *Vec4) SetVec2(in *Vec2) {
	v[0] = in.X
	v[1] = in.Y
	v[2] = 0
	v[3] = 1
}

// SetTransform overwrites the contents of this vector with the results of multiplying a provided
// vector through a provided transform matrix
//
// input - The vector to be transformed
//
// m - The transform matrix to be used in the transform
func (v *Vec4) SetTransform(input *Vec4, m *Mat4x4) {
	v[0] = m[0][0]*input[0] + m[1][0]*input[1] + m[2][0]*input[2] + m[3][0]*input[3]
	v[1] = m[0][1]*input[0] + m[1][1]*input[1] + m[2][1]*input[2] + m[3][1]*input[3]
	v[2] = m[0][2]*input[0] + m[1][2]*input[1] + m[2][2]*input[2] + m[3][2]*input[3]
	v[3] = m[0][3]*input[0] + m[1][3]*input[1] + m[2][3]*input[2] + m[3][3]*input[3]
}

// SetTransformHomogenous overwrites the contents of this vector with the results of multiplying
// a provided vector through a provided transform matrix. A homogenous divide will be performed
// on the results of the transform.
//
// input - The vector to be transformed
//
// m - The transform matrix to be used in the transform
func (v *Vec4) SetTransformHomogenous(input *Vec4, m *Mat4x4) {
	w := m[0][3]*input[0] + m[1][3]*input[1] + m[2][3]*input[2] + m[3][3]*input[3]

	v[0] = (m[0][0]*input[0] + m[1][0]*input[1] + m[2][0]*input[2] + m[3][0]*input[3]) / w
	v[1] = (m[0][1]*input[0] + m[1][1]*input[1] + m[2][1]*input[2] + m[3][1]*input[3]) / w
	v[2] = (m[0][2]*input[0] + m[1][2]*input[1] + m[2][2]*input[2] + m[3][2]*input[3]) / w
	v[3] = 1.0
}

// SetRotateWithQuaternion overwrites the contents of this vector with the results of rotating
// a provided vector through a provided quaternion.
//
// input - The vector to be rotated
//
// q - The quaternion to be used in the rotation
func (v *Vec4) SetRotateWithQuaternion(input *Vec4, q *Quaternion) {
	quatVector := Vec3{q.X, q.Y, q.Z}
	vec3 := Vec3{input[0], input[1], input[2]}
	var uv Vec3

	uv.SetCrossProduct(&quatVector, &vec3)

	var uuv Vec3
	uuv.SetCrossProduct(&quatVector, &uv)

	v[0] = input[0] + ((uv.X*q.W)+uuv.X)*2.0
	v[1] = input[1] + ((uv.Y*q.W)+uuv.Y)*2.0
	v[2] = input[2] + ((uv.Z*q.W)+uuv.Z)*2.0
	v[3] = input[3]
}

// SetRotate overwrites the contents of this vector with the results of rotating a
// provided vector around a provided axis by a provided angle
//
// input - The vector to be rotated
//
// angleRad - The amount to rotate the vector, in radians
//
// axis - The axis around which to rotate the vector
func (v *Vec4) SetRotate(input *Vec4, angleRad float64, normal *Vec3) {
	cos := float32(math.Cos(angleRad))
	sin := float32(math.Sin(angleRad))
	inverseCos := 1 - cos

	var unitAxis Vec3
	unitAxis.SetVec3(normal)
	unitAxis.Normalize()

	v[0] = (inverseCos*unitAxis.X*unitAxis.X+cos)*input[0] +
		(inverseCos*unitAxis.X*unitAxis.Y+unitAxis.Z*sin)*input[1] +
		(inverseCos*unitAxis.X*unitAxis.Z-unitAxis.Y*sin)*input[2]
	v[1] = (inverseCos*unitAxis.X*unitAxis.Y-unitAxis.Z*sin)*input[0] +
		(inverseCos*unitAxis.Y*unitAxis.Y+cos)*input[1] +
		(inverseCos*unitAxis.Y*unitAxis.Z+unitAxis.X*sin)*input[2]
	v[2] = (inverseCos*unitAxis.X*unitAxis.Z+unitAxis.Y*sin)*input[0] +
		(inverseCos*unitAxis.Y*unitAxis.Z-unitAxis.X*sin)*input[1] +
		(inverseCos*unitAxis.Z*unitAxis.Z+cos)*input[2]
	v[3] = input[3]
}

// SetRotateX overwrites the contents of this vector with the results of rotating a
// provided vector around the x axis by a provided angle
//
// input - The vector to be rotated
//
// angleRad - The amount to rotate the vector, in radians
func (v *Vec4) SetRotateX(input *Vec4, angleRad float64) {
	cos := float32(math.Cos(angleRad))
	sin := float32(math.Sin(angleRad))

	y := input[1]*cos - input[2]*sin
	z := input[1]*sin + input[2]*cos

	v[0] = input[0]
	v[1] = y
	v[2] = z
	v[3] = input[3]
}

// SetRotateY overwrites the contents of this vector with the results of rotating a
// provided vector around the y axis by a provided angle
//
// input - The vector to be rotated
//
// angleRad - The amount to rotate the vector, in radians
func (v *Vec4) SetRotateY(input *Vec4, angleRad float64) {
	cos := float32(math.Cos(angleRad))
	sin := float32(math.Sin(angleRad))

	x := input[0]*cos + input[2]*sin
	z := -input[0]*sin + input[2]*cos

	v[0] = x
	v[1] = input[1]
	v[2] = z
	v[3] = input[3]
}

// SetRotateZ overwrites the contents of this vector with the results of rotating a
// provided vector around the z axis by a provided angle
//
// input - The vector to be rotated
//
// angleRad - The amount to rotate the vector, in radians
func (v *Vec4) SetRotateZ(input *Vec4, angleRad float64) {
	cos := float32(math.Cos(angleRad))
	sin := float32(math.Sin(angleRad))

	x := input[0]*cos - input[1]*sin
	y := input[0]*sin + input[1]*cos

	v[0] = x
	v[1] = y
	v[2] = input[2]
	v[3] = input[3]
}

// SetNormalizeVec4 overwrites the contents of this vector with the results of normalizing
// a provided 4-element vector
//
// input - The vector to be normalized
func (v *Vec4) SetNormalizeVec4(input *Vec4) {
	vecLen := input.Len()

	if abs(vecLen-1) < 0.0001 {
		v[0] = input[0]
		v[1] = input[1]
		v[2] = input[2]
		v[3] = input[3]
		return
	}

	inverse := 1.0 / vecLen
	v[0] = input[0] * inverse
	v[1] = input[1] * inverse
	v[2] = input[2] * inverse
	v[3] = input[3] * inverse
}

// SetAddVec4 overwrites the contents of this vector with the results of adding two provided vectors
//
// lhs - The left operand of the cross product operation
//
// rhs - The right operand of the cross product operation
func (v *Vec4) SetAddVec4(lhs *Vec4, rhs *Vec4) {
	v[0] = lhs[0] + rhs[0]
	v[1] = lhs[1] + rhs[1]
	v[2] = lhs[2] + rhs[2]
	v[3] = lhs[3] + rhs[3]
}

// SetSubtractVec4 overwrites the contents of this vector with the results of subtracting two
// provided vectors
//
// lhs - The left operand of the cross product operation
//
// rhs - The right operand of the cross product operation
func (v *Vec4) SetSubtractVec4(lhs *Vec4, rhs *Vec4) {
	v[0] = lhs[0] - rhs[0]
	v[1] = lhs[1] - rhs[1]
	v[2] = lhs[2] - rhs[2]
	v[3] = lhs[3] - rhs[3]
}

// SetScale overwrites the contents of this vector with the provided vector
// multiplied by the provided scalar
//
// input - The vector to scale
//
// scale - The scalar t multiply with the vector
func (v *Vec4) SetScale(input *Vec4, scale float32) {
	v[0] = input[0] * scale
	v[1] = input[1] * scale
	v[2] = input[2] * scale
	v[3] = input[3] * scale
}

// SetLerp overwrites the contents of this vector with the results of linear interpolation
// between two provided vectors
//
// lhs - The source vector in the interpolation operation
//
// rhs - The target vector in the interpolation operation
//
// delta - A value between 0 and 1 indicating how far to interpolate between the two vectors
func (v *Vec4) SetLerp(lhs *Vec4, rhs *Vec4, delta float32) {
	v[0] = lhs[0]*(1-delta) + rhs[0]*delta
	v[1] = lhs[1]*(1-delta) + rhs[1]*delta
	v[2] = lhs[2]*(1-delta) + rhs[2]*delta
	v[3] = lhs[3]*(1-delta) + rhs[3]*delta
}

// Normalize converts this vector into a unit vector
func (v *Vec4) Normalize() {
	vecLen := v.Len()

	if abs(vecLen-1) < 0.0001 {
		return
	}

	inverse := 1.0 / vecLen
	v[0] *= inverse
	v[1] *= inverse
	v[2] *= inverse
	v[3] *= inverse

}

// Len calculates the length of this vector
func (v *Vec4) Len() float32 {
	sqr := float64(v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3])
	return float32(math.Sqrt(sqr))
}

// LenSqr calculates the length-squared of this vector. It is more performant than Len,
// owing to not requiring a math.Sqrt call, and may be preferable in cases when only the relative
// length of two vectors is required, or when comparing the length to 0 or 1
func (v *Vec4) LenSqr() float32 {
	return v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3]
}

// AddVec4 adds another provided vector to this one and updates this vector with the results
//
// other - The right operand in the add operation
func (v *Vec4) AddVec4(other *Vec4) {
	v[0] += other[0]
	v[1] += other[1]
	v[2] += other[2]
	v[3] += other[3]
}

// SubtractVec4 subtracts another provided vector from this one and updates this vector with the results
//
// other - The right operand in the subtract operation
func (v *Vec4) SubtractVec4(other *Vec4) {
	v[0] -= other[0]
	v[1] -= other[1]
	v[2] -= other[2]
	v[3] -= other[3]
}

// DotProduct calculates and returns the dot product of this vector and another provided vector
//
// other - The right operand in the dot product operation
func (v *Vec4) DotProduct(other *Vec4) float32 {
	return v[0]*other[0] + v[1]*other[1] + v[2]*other[2] + v[3]*other[3]
}

// Scale multiplies this vector by the provided scalar factor
//
// scale - The scalar to multiply this vector by
func (v *Vec4) Scale(scale float32) {
	v[0] *= scale
	v[1] *= scale
	v[2] *= scale
	v[3] *= scale
}

// Equal returns true if every entry in this vector is equal to every entry in the provided vector
//
// other - The vector to compare this vector to
//
// epsilon - The epsilon value to use in floating point comparisons. This much floating point
// drift is permitted before the method returns false. 0.0001 is a common epsilon value
func (v *Vec4) Equal(other *Vec4, epsilon float32) bool {
	if abs(v[0]-other[0]) > epsilon {
		return false
	}

	if abs(v[1]-other[1]) > epsilon {
		return false
	}

	if abs(v[2]-other[2]) > epsilon {
		return false
	}

	if abs(v[3]-other[3]) > epsilon {
		return false
	}

	return true
}

// Transform multiplies this vector through the provided transform matrix
// and updates the vector with the result
//
// m - The 4x4 matrix to transform this matrix through
func (v *Vec4) Transform(m *Mat4x4) {
	TransformVec(v, m, v)
}

// TransformHomogenous performs a transform with a provided 4x4 matrix. After the transform,
// a homogenous divide is performed on the results. This vector is updated with the output.
//
// m - The 4x4 matrix to transform this matrix through
func (v *Vec4) TransformHomogenous(m *Mat4x4) {
	x := m[0][0]*v[0] + m[1][0]*v[1] + m[2][0]*v[2] + m[3][0]*v[3]
	y := m[0][1]*v[0] + m[1][1]*v[1] + m[2][1]*v[2] + m[3][1]*v[3]
	z := m[0][2]*v[0] + m[1][2]*v[1] + m[2][2]*v[2] + m[3][2]*v[3]
	w := m[0][3]*v[0] + m[1][3]*v[1] + m[2][3]*v[2] + m[3][3]*v[3]

	v[0] = x / w
	v[1] = y / w
	v[2] = z / w
	v[3] = 1.0
}

// Lerp interpolates between this vector and another provided vector, updating this
// vector to the results of the interpolation
//
// other - The target vector in the interpolation operation
//
// delta - A value between 0 and 1 indicating how far to interpolate between the two vectors
func (v *Vec4) Lerp(other *Vec4, delta float32) {
	v[0] = v[0]*(1-delta) + other[0]*delta
	v[1] = v[1]*(1-delta) + other[1]*delta
	v[2] = v[2]*(1-delta) + other[2]*delta
	v[3] = v[3]*(1-delta) + other[3]*delta
}

// Rotate rotates this vector around the provided axis by the provided angle
//
// angleRad - The amount to rotate this vector, in radians
//
// axis - The axis to rotate this vector around
func (v *Vec4) Rotate(angleRad float64, normal *Vec3) {
	cos := float32(math.Cos(angleRad))
	sin := float32(math.Sin(angleRad))
	inverseCos := 1 - cos

	var unitAxis Vec3
	unitAxis.SetVec3(normal)
	unitAxis.Normalize()

	x := (inverseCos*unitAxis.X*unitAxis.X+cos)*v[0] +
		(inverseCos*unitAxis.X*unitAxis.Y+unitAxis.Z*sin)*v[1] +
		(inverseCos*unitAxis.X*unitAxis.Z-unitAxis.Y*sin)*v[2]

	y := (inverseCos*unitAxis.X*unitAxis.Y-unitAxis.Z*sin)*v[0] +
		(inverseCos*unitAxis.Y*unitAxis.Y+cos)*v[1] +
		(inverseCos*unitAxis.Y*unitAxis.Z+unitAxis.X*sin)*v[2]

	z := (inverseCos*unitAxis.X*unitAxis.Z+unitAxis.Y*sin)*v[0] +
		(inverseCos*unitAxis.Y*unitAxis.Z-unitAxis.X*sin)*v[1] +
		(inverseCos*unitAxis.Z*unitAxis.Z+cos)*v[2]

	v[0] = x
	v[1] = y
	v[2] = z
}

// RotateWithQuaternion rotates this vector with the provided quaternion
//
// q - The quaternion to rotate this vector with
func (v *Vec4) RotateWithQuaternion(q *Quaternion) {
	quatVector := Vec3{q.X, q.Y, q.Z}
	vec3 := Vec3{v[0], v[1], v[2]}

	var uv Vec3
	uv.SetCrossProduct(&quatVector, &vec3)

	var uuv Vec3
	uuv.SetCrossProduct(&quatVector, &uv)

	v[0] += ((uv.X * q.W) + uuv.X) * 2.0
	v[1] += ((uv.Y * q.W) + uuv.Y) * 2.0
	v[2] += ((uv.Z * q.W) + uuv.Z) * 2.0
}

// RotateX rotates this vector around the x axis by a provided angle
//
// angleRad - The amount to rotate the vector, in radians
func (v *Vec4) RotateX(angleRad float64) {
	cos := float32(math.Cos(angleRad))
	sin := float32(math.Sin(angleRad))

	y := v[1]*cos - v[2]*sin
	z := v[1]*sin + v[2]*cos

	v[1] = y
	v[2] = z
}

// RotateY rotates this vector around the y axis by a provided angle
//
// angleRad - The amount to rotate the vector, in radians
func (v *Vec4) RotateY(angleRad float64) {
	cos := float32(math.Cos(angleRad))
	sin := float32(math.Sin(angleRad))

	x := v[0]*cos + v[2]*sin
	z := -v[0]*sin + v[2]*cos

	v[0] = x
	v[2] = z
}

// RotateZ rotates this vector around the z axis by a provided angle
//
// angleRad - The amount to rotate the vector, in radians
func (v *Vec4) RotateZ(angleRad float64) {
	cos := float32(math.Cos(angleRad))
	sin := float32(math.Sin(angleRad))

	x := v[0]*cos - v[1]*sin
	y := v[0]*sin + v[1]*cos

	v[0] = x
	v[1] = y
}
