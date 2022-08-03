package math

import "math"

type Vec4[T FloatingPoint] struct {
	X T
	Y T
	Z T
	W T
}

func (v *Vec4[T]) CopyFrom(in *Vec4[T]) *Vec4[T] {
	v.X = in.X
	v.Y = in.Y
	v.Z = in.Z
	v.W = in.W

	return v
}

func (v *Vec4[T]) Normalize() *Vec4[T] {
	inverse := 1.0 / v.Len()
	v.X *= inverse
	v.Y *= inverse
	v.Z *= inverse
	v.W *= inverse

	return v
}

func (v *Vec4[T]) Len() T {
	sqr := float64(v.X*v.X + v.Y*v.Y + v.Z*v.Z + v.W*v.W)
	return T(math.Sqrt(sqr))
}

func (v *Vec4[T]) LenSqr() T {
	return v.X*v.X + v.Y*v.Y + v.Z*v.Z + v.W*v.W
}

func (v *Vec4[T]) AddVec4(other *Vec4[T]) *Vec4[T] {
	v.X += other.X
	v.Y += other.Y
	v.Z += other.Z
	v.W += other.W

	return v
}

func (v *Vec4[T]) SubtractVec4(other *Vec4[T]) *Vec4[T] {
	v.X -= other.X
	v.Y -= other.Y
	v.Z -= other.Z
	v.W -= other.W

	return v
}

func (v *Vec4[T]) DotProduct(other *Vec4[T]) T {
	return v.X*other.X + v.Y*other.Y + v.Z*other.Z + v.W*other.W
}

func (v *Vec4[T]) Scale(scale T) *Vec4[T] {
	v.X *= scale
	v.Y *= scale
	v.Z *= scale
	v.W *= scale

	return v
}

func (v *Vec4[T]) Equal(other *Vec4[T], epsilon T) bool {
	if abs[T](v.X-other.X) > epsilon {
		return false
	}

	if abs[T](v.Y-other.Y) > epsilon {
		return false
	}

	if abs[T](v.Z-other.Z) > epsilon {
		return false
	}

	if abs[T](v.W-other.W) > epsilon {
		return false
	}

	return true
}
