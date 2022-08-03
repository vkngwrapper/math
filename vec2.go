package math

import "math"

type Vec2[T FloatingPoint] struct {
	X T
	Y T
}

func (v *Vec2[T]) CopyFrom(in *Vec3[T]) *Vec2[T] {
	v.X = in.X
	v.Y = in.Y
	return v
}

func (v *Vec2[T]) Normalize() *Vec2[T] {
	inverse := 1.0 / v.Len()
	v.X *= inverse
	v.Y *= inverse

	return v
}

func (v *Vec2[T]) Len() T {
	sqr := float64(v.X*v.X + v.Y*v.Y)
	return T(math.Sqrt(sqr))
}

func (v *Vec2[T]) LenSqr() T {
	return v.X*v.X + v.Y*v.Y
}

func (v *Vec2[T]) AddVec3(other *Vec2[T]) *Vec2[T] {
	v.X += other.X
	v.Y += other.Y

	return v
}

func (v *Vec2[T]) SubtractVec3(other *Vec2[T]) *Vec2[T] {
	v.X -= other.X
	v.Y -= other.Y

	return v
}

func (v *Vec2[T]) DotProduct(other *Vec2[T]) T {
	return v.X*other.X + v.Y*other.Y
}

func (v *Vec2[T]) Scale(scale T) *Vec2[T] {
	v.X *= scale
	v.Y *= scale

	return v
}

func (v *Vec2[T]) Equal(other *Vec2[T], epsilon T) bool {
	if abs[T](v.X-other.X) > epsilon {
		return false
	}

	if abs[T](v.Y-other.Y) > epsilon {
		return false
	}

	return true
}
