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

func (v *Vec4[T]) Transform(m *Mat4x4[T]) *Vec4[T] {
	x := m[0][0]*v.X + m[1][0]*v.Y + m[2][0]*v.Z + m[3][0]*v.W
	y := m[0][1]*v.X + m[1][1]*v.Y + m[2][1]*v.Z + m[3][1]*v.W
	z := m[0][2]*v.X + m[1][2]*v.Y + m[2][2]*v.Z + m[3][2]*v.W
	w := m[0][3]*v.X + m[1][3]*v.Y + m[2][3]*v.Z + m[3][3]*v.W

	v.X = x
	v.Y = y
	v.Z = z
	v.W = w

	return v
}

func (v *Vec4[T]) TransformExp(m *Mat4x4Exp[T]) *Vec4[T] {
	x := m.C00*v.X + m.C10*v.Y + m.C20*v.Z + m.C30*v.W
	y := m.C01*v.X + m.C11*v.Y + m.C21*v.Z + m.C31*v.W
	z := m.C02*v.X + m.C12*v.Y + m.C22*v.Z + m.C32*v.W
	w := m.C03*v.X + m.C13*v.Y + m.C23*v.Z + m.C33*v.W

	v.X = x
	v.Y = y
	v.Z = z
	v.W = w

	return v
}
