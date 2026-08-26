package math

import "math"

// Vec2 is a two-element vector of floating point compatible values that can be used for 2d
// vector arithmetic
type Vec2[T FloatingPoint] struct {
	X T
	Y T
}

// SetVec2 overwrites the contents of this vector with the contents of a 2-element vector
//
// in - The vector to initialize from
func (v *Vec2[T]) SetVec2(in *Vec2[T]) {
	v.X = in.X
	v.Y = in.Y
}

// SetVec3 overwrites the contents of this vector with the first two elements of a 3-element vector
//
// in - The vector to initialize from
func (v *Vec2[T]) SetVec3(in *Vec3[T]) {
	v.X = in.X
	v.Y = in.Y
}

// SetVec4 overwrites the contents of this vector with the first two elements of a 4-element vector
//
// in - The vector to initialize from
func (v *Vec2[T]) SetVec4(in *Vec4[T]) {
	v.X = in.X
	v.Y = in.Y
}

// SetHomogenousVec3 overwrites the contents of this vector with the first two elements of a 3-element
// vector, divided by the third element
//
// in - The vector to initialize from
func (v *Vec2[T]) SetHomogenousVec3(in *Vec3[T]) {
	factor := T(1) / in.Z

	v.X = in.X * factor
	v.Y = in.Y * factor
}

// SetHomogenousVec4 overwrites the contents of this vector with the first two elements of a 4-element
// vector, divided by the fourth element
//
// in - The vector to initialize from
func (v *Vec2[T]) SetHomogenousVec4(in *Vec4[T]) {
	factor := T(1) / in.W

	v.X = in.X * factor
	v.Y = in.Y * factor
}

// Normalize converts this vector into a unit vector
func (v *Vec2[T]) Normalize() {
	vecLen := v.Len()

	if abs[T](vecLen-1) > 0.0001 {
		return
	}

	inverse := 1.0 / vecLen
	v.X *= inverse
	v.Y *= inverse
}

// Len calculates the length of this vector
func (v *Vec2[T]) Len() T {
	sqr := float64(v.X*v.X + v.Y*v.Y)
	return T(math.Sqrt(sqr))
}

// LenSqr calculates the length-squared of this vector. It is more performant than Len,
// owing to not requiring a math.Sqrt call, and may be preferable in cases when only the relative
// length of two vectors is required, or when comparing the length to 0 or 1
func (v *Vec2[T]) LenSqr() T {
	return v.X*v.X + v.Y*v.Y
}

// AddVec2 adds another provided vector to this one and updates this vector with the results
//
// other - The right operand in the add operation
func (v *Vec2[T]) AddVec2(other *Vec2[T]) {
	v.X += other.X
	v.Y += other.Y
}

// SubtractVec2 subtracts another provided vector from this one and updates this vector with the results
//
// other - The right operand in the subtract operation
func (v *Vec2[T]) SubtractVec2(other *Vec2[T]) {
	v.X -= other.X
	v.Y -= other.Y
}

// DotProduct calculates and returns the dot product of this vector and another provided vector
//
// other - The right operand in the dot product operation
func (v *Vec2[T]) DotProduct(other *Vec2[T]) T {
	return v.X*other.X + v.Y*other.Y
}

// Scale multiplies this vector by the provided scalar factor
//
// scale - The scalar to multiply this vector by
func (v *Vec2[T]) Scale(scale T) {
	v.X *= scale
	v.Y *= scale
}

// Equal returns true if every entry in this vector is equal to every entry in the provided vector
//
// other - The vector to compare this vector to
//
// epsilon - The epsilon value to use in floating point comparisons. This much floating point
// drift is permitted before the method returns false. 0.0001 is a common epsilon value
func (v *Vec2[T]) Equal(other *Vec2[T], epsilon T) bool {
	if abs[T](v.X-other.X) > epsilon {
		return false
	}

	if abs[T](v.Y-other.Y) > epsilon {
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
func (v *Vec2[T]) Lerp(other *Vec2[T], delta T) {
	v.X = v.X*(1-delta) + other.X*delta
	v.Y = v.Y*(1-delta) + other.Y*delta
}

// Rotate rotates this vector around the origin by the provided angle
//
// angleRad - The angle to rotate, in radians
func (v *Vec2[T]) Rotate(angleRad float64) {
	cos := T(math.Cos(angleRad))
	sin := T(math.Sin(angleRad))

	x := v.X*cos - v.Y*sin
	y := v.X*sin + v.Y*cos

	v.X = x
	v.Y = y
}

// Transform multiplies this vector through the provided 2x2 transform matrix
//
// m - The transform matrix used to transform this vector
func (v *Vec2[T]) Transform(m *Mat2x2[T]) {
	x := m[0][0]*v.X + m[1][0]*v.Y
	y := m[0][1]*v.X + m[1][1]*v.Y

	v.X = x
	v.Y = y
}

// TransformHomogenous multiplies this vector through the provided 3x3 transform matrix
// and updates this 2d vector with the results by performing a homogenous divide
//
// m - The transform matrix used to transform this vector
func (v *Vec2[T]) TransformHomogenous(m *Mat3x3[T]) {
	x := m[0][0]*v.X + m[1][0]*v.Y + m[2][0]
	y := m[0][1]*v.X + m[1][1]*v.Y + m[2][1]
	z := m[0][2]*v.X + m[1][2]*v.Y + m[2][2]

	factor := 1.0 / z
	v.X = x * factor
	v.Y = y * factor
}
