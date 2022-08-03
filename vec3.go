package math

import "math"

type Vec3[T FloatingPoint] struct {
	X T
	Y T
	Z T
}

func (v *Vec3[T]) CopyFrom(in *Vec3[T]) *Vec3[T] {
	v.X = in.X
	v.Y = in.Y
	v.Z = in.Z

	return v
}

func (v *Vec3[T]) Normalize() *Vec3[T] {
	inverse := 1.0 / v.Len()
	v.X *= inverse
	v.Y *= inverse
	v.Z *= inverse

	return v
}

func (v *Vec3[T]) Len() T {
	sqr := float64(v.X*v.X + v.Y*v.Y + v.Z*v.Z)
	return T(math.Sqrt(sqr))
}

func (v *Vec3[T]) LenSqr() T {
	return v.X*v.X + v.Y*v.Y + v.Z*v.Z
}

func (v *Vec3[T]) AddVec3(other *Vec3[T]) *Vec3[T] {
	v.X += other.X
	v.Y += other.Y
	v.Z += other.Z

	return v
}

func (v *Vec3[T]) SubtractVec3(other *Vec3[T]) *Vec3[T] {
	v.X -= other.X
	v.Y -= other.Y
	v.Z -= other.Z

	return v
}

func (v *Vec3[T]) CrossProduct(other *Vec3[T]) *Vec3[T] {
	x := v.Y*other.Z - v.Z*other.Y
	y := v.Z*other.X - v.X*other.Z
	z := v.X*other.Y - v.Y*other.X

	v.X = x
	v.Y = y
	v.Z = z

	return v
}

func (v *Vec3[T]) DotProduct(other *Vec3[T]) T {
	return v.X*other.X + v.Y*other.Y + v.Z*other.Z
}

func (v *Vec3[T]) Scale(scale T) *Vec3[T] {
	v.X *= scale
	v.Y *= scale
	v.Z *= scale

	return v
}
