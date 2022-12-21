package math

import (
	"math"
)

// Vec3 is a three-element vector of floating point compatible values that can be used for 3d
// vector arithmetic and transforms
type Vec3[T FloatingPoint] struct {
	X T
	Y T
	Z T
}

// SetVec3 overwrites the contents of this vector with the contents of a 3-element vector
//
// in - The vector to initialize from
func (v *Vec3[T]) SetVec3(in *Vec3[T]) *Vec3[T] {
	v.X = in.X
	v.Y = in.Y
	v.Z = in.Z

	return v
}

// SetVec2 overwrites the first two elements of this vector with the first two elements of a 2-element vector.
// The third element of this vector is set to 0.
//
// in - The vector to initialize from
func (v *Vec3[T]) SetVec2(in *Vec2[T]) *Vec3[T] {
	v.X = in.X
	v.Y = in.Y
	v.Z = 0

	return v
}

// SetVec4 overwrites the contents of this vector with the first three elements of a 4-element vector
//
// in - The vector to initialize from
func (v *Vec3[T]) SetVec4(in *Vec4[T]) *Vec3[T] {
	v.X = in.X
	v.Y = in.Y
	v.Z = in.Z

	return v
}

// SetHomogenousVec4 overwrites the contents of this vector with the first three elements of a 4-element
// vector, divided by the fourth element
//
// in - The vector to initialize from
func (v *Vec3[T]) SetHomogenousVec4(in *Vec4[T]) *Vec3[T] {
	factor := T(1) / in.W

	v.X = in.X * factor
	v.Y = in.Y * factor
	v.Z = in.Z * factor

	return v
}

// SetCrossProduct overwrites the contents of this vector with the cross product of two provided
// vectors
//
// lhs - The left operand of the cross product operation
//
// rhs - The right operand of the cross product operation
func (v *Vec3[T]) SetCrossProduct(lhs, rhs *Vec3[T]) *Vec3[T] {
	v.X = lhs.Y*rhs.Z - lhs.Z*rhs.Y
	v.Y = lhs.Z*rhs.X - lhs.X*rhs.Z
	v.Z = lhs.X*rhs.Y - lhs.Y*rhs.X

	return v
}

// SetSubtractVec3 overwrites the contents of this vector with the results of subtracting two
// provided vectors
//
// lhs - The left operand of the cross product operation
//
// rhs - The right operand of the cross product operation
func (v *Vec3[T]) SetSubtractVec3(lhs, rhs *Vec3[T]) *Vec3[T] {
	v.X = lhs.X - rhs.X
	v.Y = lhs.Y - rhs.Y
	v.Z = lhs.Z - rhs.Z

	return v
}

// SetAddVec3 overwrites the contents of this vector with the results of adding two provided vectors
//
// lhs - The left operand of the cross product operation
//
// rhs - The right operand of the cross product operation
func (v *Vec3[T]) SetAddVec3(lhs, rhs *Vec3[T]) *Vec3[T] {
	v.X = lhs.X + rhs.X
	v.Y = lhs.Y + rhs.Y
	v.Z = lhs.Z + rhs.Z

	return v
}

// SetTransform overwrites the contents of this vector with the results of multiplying a provided
// vector through a provided 3x3 transform matrix
//
// input - The vector to be transformed
//
// m - The transform matrix to be used in the transform
func (v *Vec3[T]) SetTransform(input *Vec3[T], m *Mat3x3[T]) *Vec3[T] {
	x := m[0][0]*input.X + m[1][0]*input.Y + m[2][0]*input.Z
	y := m[0][1]*input.X + m[1][1]*input.Y + m[2][1]*input.Z
	z := m[0][2]*input.X + m[1][2]*input.Y + m[2][2]*input.Z

	v.X = x
	v.Y = y
	v.Z = z

	return v
}

// SetTransformHomogenous overwrites the contents of this vector with the results of multiplying
// a provided 3-element vector through a provided 4x4 transform matrix. During the multiplication,
// the value 1 will be used as a fourth vector element, and when overwriting the vector contents,
// the 3 assigned elements will be divided by the 4th element output from the transform.
//
// input - The vector to be transformed
//
// m - The transform matrix to be used in the transform
func (v *Vec3[T]) SetTransformHomogenous(input *Vec3[T], m *Mat4x4[T]) *Vec3[T] {
	x := m[0][0]*input.X + m[1][0]*input.Y + m[2][0]*input.Z + m[3][0]
	y := m[0][1]*input.X + m[1][1]*input.Y + m[2][1]*input.Z + m[3][1]
	z := m[0][2]*input.X + m[1][2]*input.Y + m[2][2]*input.Z + m[3][2]
	w := m[0][3]*input.X + m[1][3]*input.Y + m[2][3]*input.Z + m[3][3]

	factor := 1.0 / w
	v.X = x * factor
	v.Y = y * factor
	v.Z = z * factor

	return v
}

// SetRotateWithQuaternion overwrites the contents of this vector with the results of rotating
// a provided 3-element vector through a provided quaternion.
//
// input - The vector to be rotated
//
// q - The quaternion to be used in the rotation
func (v *Vec3[T]) SetRotateWithQuaternion(input *Vec3[T], q *Quaternion[T]) *Vec3[T] {
	quatVector := Vec3[T]{q.X, q.Y, q.Z}

	var uv Vec3[T]
	uv.SetCrossProduct(&quatVector, input)

	var uuv Vec3[T]
	uuv.SetCrossProduct(&quatVector, &uv)

	v.X = input.X + ((uv.X*q.W)+uuv.X)*2.0
	v.Y = input.Y + ((uv.Y*q.W)+uuv.Y)*2.0
	v.Z = input.Z + ((uv.Z*q.W)+uuv.Z)*2.0

	return v
}

// SetRotate overwrites the contents of this vector with the results of rotating a
// provided 3-element vector around a provided axis by a provided angle
//
// input - The vector to be rotated
//
// angleRad - The amount to rotate the vector, in radians
//
// axis - The axis around which to rotate the vector
func (v *Vec3[T]) SetRotate(input *Vec3[T], angleRad float64, axis *Vec3[T]) *Vec3[T] {
	cos := T(math.Cos(angleRad))
	sin := T(math.Sin(angleRad))
	inverseCos := T(1) - cos

	var unitAxis Vec3[T]
	unitAxis.SetVec3(axis).Normalize()

	v.X = (inverseCos*unitAxis.X*unitAxis.X+cos)*input.X +
		(inverseCos*unitAxis.X*unitAxis.Y+unitAxis.Z*sin)*input.Y +
		(inverseCos*unitAxis.X*unitAxis.Z-unitAxis.Y*sin)*input.Z
	v.Y = (inverseCos*unitAxis.X*unitAxis.Y-unitAxis.Z*sin)*input.X +
		(inverseCos*unitAxis.Y*unitAxis.Y+cos)*input.Y +
		(inverseCos*unitAxis.Y*unitAxis.Z+unitAxis.X*sin)*input.Z
	v.Z = (inverseCos*unitAxis.X*unitAxis.Z+unitAxis.Y*sin)*input.X +
		(inverseCos*unitAxis.Y*unitAxis.Z-unitAxis.X*sin)*input.Y +
		(inverseCos*unitAxis.Z*unitAxis.Z+cos)*input.Z

	return v
}

// SetRotateX overwrites the contents of this vector with the results of rotating a
// provided 3-element vector around the x axis by a provided angle
//
// input - The vector to be rotated
//
// angleRad - The amount to rotate the vector, in radians
func (v *Vec3[T]) SetRotateX(input *Vec3[T], angleRad float64) *Vec3[T] {
	cos := T(math.Cos(angleRad))
	sin := T(math.Sin(angleRad))

	v.X = input.X
	v.Y = input.Y*cos - input.Z*sin
	v.Z = input.Y*sin + input.Z*cos

	return v
}

// SetRotateY overwrites the contents of this vector with the results of rotating a
// provided 3-element vector around the y axis by a provided angle
//
// input - The vector to be rotated
//
// angleRad - The amount to rotate the vector, in radians
func (v *Vec3[T]) SetRotateY(input *Vec3[T], angleRad float64) *Vec3[T] {
	cos := T(math.Cos(angleRad))
	sin := T(math.Sin(angleRad))

	v.X = input.X*cos + input.Z*sin
	v.Y = input.Y
	v.Z = -input.X*sin + input.Z*cos

	return v
}

// SetRotateZ overwrites the contents of this vector with the results of rotating a
// provided 3-element vector around the z axis by a provided angle
//
// input - The vector to be rotated
//
// angleRad - The amount to rotate the vector, in radians
func (v *Vec3[T]) SetRotateZ(input *Vec3[T], angleRad float64) *Vec3[T] {
	cos := T(math.Cos(angleRad))
	sin := T(math.Sin(angleRad))

	v.X = input.X*cos - input.Y*sin
	v.Y = input.X*sin + input.Y*cos
	v.Z = input.Z

	return v
}

// SetNormalizeVec3 overwrites the contents of this vector with the results of normalizing
// a provided 3-element vector
//
// input - The vector to be normalized
func (v *Vec3[T]) SetNormalizeVec3(input *Vec3[T]) *Vec3[T] {
	vecLen := input.Len()

	if abs[T](vecLen-1) < 0.0001 {
		v.X = input.X
		v.Y = input.Y
		v.Z = input.Z
		return v
	}

	inverse := 1.0 / vecLen
	v.X = input.X * inverse
	v.Y = input.Y * inverse
	v.Z = input.Z * inverse

	return v
}

// SetLerp overwrites the contents of this vector with the results of linear interpolation
// between two provided 3-element vectors
//
// lhs - The source vector in the interpolation operation
//
// rhs - The target vector in the interpolation operation
//
// delta - A value between 0 and 1 indicating how far to interpolate between the two vectors
func (v *Vec3[T]) SetLerp(lhs, rhs *Vec3[T], delta T) *Vec3[T] {
	v.X = lhs.X*(1-delta) + rhs.X*delta
	v.Y = lhs.Y*(1-delta) + rhs.Y*delta
	v.Z = lhs.Z*(1-delta) + rhs.Z*delta

	return v
}

// Normalize converts this vector into a unit vector
func (v *Vec3[T]) Normalize() *Vec3[T] {
	vecLenSqr := v.LenSqr()

	if abs[T](vecLenSqr-1) < 0.0001 {
		return v
	}

	vecLen := T(math.Sqrt(float64(vecLenSqr)))
	inverse := 1.0 / vecLen
	v.X *= inverse
	v.Y *= inverse
	v.Z *= inverse

	return v
}

// Len calculates the length of this vector
func (v *Vec3[T]) Len() T {
	sqr := float64(v.X*v.X + v.Y*v.Y + v.Z*v.Z)
	return T(math.Sqrt(sqr))
}

// LenSqr calculates the length-squared of this vector. It is more performant than Len,
// owing to not requiring a math.Sqrt call, and may be preferable in cases when only the relative
// length of two vectors is required, or when comparing the length to 0 or 1
func (v *Vec3[T]) LenSqr() T {
	return v.X*v.X + v.Y*v.Y + v.Z*v.Z
}

// AddVec3 adds another provided vector to this one and updates this vector with the results
//
// other - The right operand in the add operation
func (v *Vec3[T]) AddVec3(other *Vec3[T]) *Vec3[T] {
	v.X += other.X
	v.Y += other.Y
	v.Z += other.Z

	return v
}

// SubtractVec3 subtracts another provided vector from this one and updates this vector with the results
//
// other - The right operand in the subtract operation
func (v *Vec3[T]) SubtractVec3(other *Vec3[T]) *Vec3[T] {
	v.X -= other.X
	v.Y -= other.Y
	v.Z -= other.Z

	return v
}

// CrossProduct calculates the cross product of this vector and another provided vector and updates this
// vector to the calculated cross product
//
// other - The right operand of the cross product operation
func (v *Vec3[T]) CrossProduct(other *Vec3[T]) *Vec3[T] {
	x := v.Y*other.Z - v.Z*other.Y
	y := v.Z*other.X - v.X*other.Z
	z := v.X*other.Y - v.Y*other.X

	v.X = x
	v.Y = y
	v.Z = z

	return v
}

// DotProduct calculates and returns the dot product of this vector and another provided vector
//
// other - The right operand in the dot product operation
func (v *Vec3[T]) DotProduct(other *Vec3[T]) T {
	return v.X*other.X + v.Y*other.Y + v.Z*other.Z
}

// SetScale overwrites the contents of this vector with the provided vector
// multiplied by the provided scalar
//
// input - The vector to scale
//
// scale - The scalar t multiply with the vector
func (v *Vec3[T]) SetScale(input *Vec3[T], scale T) *Vec3[T] {
	v.X = input.X * scale
	v.Y = input.Y * scale
	v.Z = input.Z * scale

	return v
}

// Scale multiplies this vector by the provided scalar factor
//
// scale - The scalar to multiply this vector by
func (v *Vec3[T]) Scale(scale T) *Vec3[T] {
	v.X *= scale
	v.Y *= scale
	v.Z *= scale

	return v
}

// Equal returns true if every entry in this vector is equal to every entry in the provided vector
//
// other - The vector to compare this vector to
//
// epsilon - The epsilon value to use in floating point comparisons. This much floating point
// drift is permitted before the method returns false. 0.0001 is a common epsilon value
func (v *Vec3[T]) Equal(other *Vec3[T], epsilon T) bool {
	if abs[T](v.X-other.X) > epsilon {
		return false
	}

	if abs[T](v.Y-other.Y) > epsilon {
		return false
	}

	if abs[T](v.Z-other.Z) > epsilon {
		return false
	}

	return true
}

// Lerp interpolates between this vector and another provided vector, updating this
// vector to the results of the interpolation
//
// other - The target vector in the interpolation operation
//
// delta - A value between 0 and 1 indicating how far to interpolate between the two vectors
func (v *Vec3[T]) Lerp(other *Vec3[T], delta T) *Vec3[T] {
	v.X = v.X*(1-delta) + other.X*delta
	v.Y = v.Y*(1-delta) + other.Y*delta
	v.Z = v.Z*(1-delta) + other.Z*delta

	return v
}

// SetClosestPointOnLine overwrites the contents of this vector with the point on the line
// that includes both endpoint1 and endpoint2 that falls closest to the provided point
//
// point - The point to locate the closest position to
//
// endpoint1 - The first point of the line the result will fall along
//
// endpoint2 - The second point of the line the result will fall along
func (v *Vec3[T]) SetClosestPointOnLine(point, endpoint1, endpoint2 *Vec3[T]) *Vec3[T] {

	var unitLine Vec3[T]
	unitLine.SetSubtractVec3(endpoint2, endpoint1)

	distanceSquared := unitLine.LenSqr()

	unitLine.Normalize()

	var lineToPoint Vec3[T]
	lineToPoint.SetSubtractVec3(point, endpoint1)

	distance := lineToPoint.DotProduct(&unitLine)

	if distance <= 0 {
		return endpoint1
	} else if distance*distance >= distanceSquared {
		return endpoint2
	}

	v.X = endpoint1.X + unitLine.X*distance
	v.Y = endpoint1.Y + unitLine.Y*distance
	v.Z = endpoint1.Z + unitLine.Z*distance

	return v
}

// Orthonormalize adjusts this vector to be a unit vector that lies orthogonal to the provided vector
//
// other - The vector that this vector will lie orthogonal to after the operation
func (v *Vec3[T]) Orthonormalize(other *Vec3[T]) *Vec3[T] {
	dotProduct := v.DotProduct(other)

	v.X = v.X - other.X*dotProduct
	v.Y = v.Y - other.Y*dotProduct
	v.Z = v.Z - other.Z*dotProduct

	return v.Normalize()
}

// Rotate rotates this vector around the provided axis by the provided angle
//
// angleRad - The amount to rotate this vector, in radians
//
// axis - The axis to rotate this vector around
func (v *Vec3[T]) Rotate(angleRad float64, axis *Vec3[T]) *Vec3[T] {
	cos := T(math.Cos(angleRad))
	sin := T(math.Sin(angleRad))
	inverseCos := T(1) - cos

	var unitAxis Vec3[T]
	unitAxis.SetVec3(axis).Normalize()
	dotProduct := v.DotProduct(&unitAxis)

	x := cos*v.X + (unitAxis.Y*v.Z-unitAxis.Z*v.Y)*sin + unitAxis.X*dotProduct*inverseCos
	y := cos*v.Y + (unitAxis.Z*v.X-unitAxis.X*v.Z)*sin + unitAxis.Y*dotProduct*inverseCos
	z := cos*v.Z + (unitAxis.X*v.Y-unitAxis.Y*v.X)*sin + unitAxis.Z*dotProduct*inverseCos

	v.X = x
	v.Y = y
	v.Z = z

	return v
}

// RotateWithQuaternion rotates this vector with the provided quaternion
//
// q - The quaternion to rotate this vector with
func (v *Vec3[T]) RotateWithQuaternion(q *Quaternion[T]) *Vec3[T] {
	quatVector := Vec3[T]{q.X, q.Y, q.Z}

	var uv Vec3[T]
	uv.SetCrossProduct(&quatVector, v)

	var uuv Vec3[T]
	uuv.SetCrossProduct(&quatVector, &uv)

	v.X += ((uv.X * q.W) + uuv.X) * 2.0
	v.Y += ((uv.Y * q.W) + uuv.Y) * 2.0
	v.Z += ((uv.Z * q.W) + uuv.Z) * 2.0

	return v
}

// RotateX rotates this vector around the x axis by a provided angle
//
// angleRad - The amount to rotate the vector, in radians
func (v *Vec3[T]) RotateX(angleRad float64) *Vec3[T] {
	cos := T(math.Cos(angleRad))
	sin := T(math.Sin(angleRad))

	y := v.Y*cos - v.Z*sin
	z := v.Y*sin + v.Z*cos

	v.Y = y
	v.Z = z

	return v
}

// RotateY rotates this vector around the y axis by a provided angle
//
// angleRad - The amount to rotate the vector, in radians
func (v *Vec3[T]) RotateY(angleRad float64) *Vec3[T] {
	cos := T(math.Cos(angleRad))
	sin := T(math.Sin(angleRad))

	x := v.X*cos + v.Z*sin
	z := -v.X*sin + v.Z*cos

	v.X = x
	v.Z = z

	return v
}

// RotateZ rotates this vector around the z axis by a provided angle
//
// angleRad - The amount to rotate the vector, in radians
func (v *Vec3[T]) RotateZ(angleRad float64) *Vec3[T] {
	cos := T(math.Cos(angleRad))
	sin := T(math.Sin(angleRad))

	x := v.X*cos - v.Y*sin
	y := v.X*sin + v.Y*cos

	v.X = x
	v.Y = y

	return v
}

// SetTriangleNormal overwrites the contents of this vector with the unit normal of a triangle represented
// by the provided points. This method assumes a counter-clockwise winding order.
//
// point1 - The first point of the triangle
//
// point2 - The second point of the triangle
//
// point3 - The third point of the triangle
func (v *Vec3[T]) SetTriangleNormal(point1, point2, point3 *Vec3[T]) *Vec3[T] {
	var edge1, edge2 Vec3[T]

	edge2.SetSubtractVec3(point1, point3)
	normal := edge1.SetSubtractVec3(point1, point2).CrossProduct(&edge2).Normalize()

	v.X = normal.X
	v.Y = normal.Y
	v.Z = normal.Z

	return v
}

// Transform multiplies this vector through the provided transform matrix
// and updates the vector with the result
//
// m - The 3x3 matrix to transform this matrix through
func (v *Vec3[T]) Transform(m *Mat3x3[T]) *Vec3[T] {
	x := m[0][0]*v.X + m[1][0]*v.Y + m[2][0]*v.Z
	y := m[0][1]*v.X + m[1][1]*v.Y + m[2][1]*v.Z
	z := m[0][2]*v.X + m[1][2]*v.Y + m[2][2]*v.Z

	v.X = x
	v.Y = y
	v.Z = z

	return v
}

// TransformHomogenous performs a transform with a provided 4x4 matrix. During the transform,
// the value 1 is provided as the 4th vector element and a homogenous divide is performed on
// the results. This vector is updated with the output.
//
// m - The 4x4 matrix to transform this matrix through
func (v *Vec3[T]) TransformHomogenous(m *Mat4x4[T]) *Vec3[T] {
	x := m[0][0]*v.X + m[1][0]*v.Y + m[2][0]*v.Z + m[3][0]
	y := m[0][1]*v.X + m[1][1]*v.Y + m[2][1]*v.Z + m[3][1]
	z := m[0][2]*v.X + m[1][2]*v.Y + m[2][2]*v.Z + m[3][2]
	w := m[0][3]*v.X + m[1][3]*v.Y + m[2][3]*v.Z + m[3][3]

	factor := 1.0 / w
	v.X = x * factor
	v.Y = y * factor
	v.Z = z * factor

	return v
}
