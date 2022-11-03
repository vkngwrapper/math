package math

import "math"

type Vec2[T FloatingPoint] struct {
	X T
	Y T
}

func (v *Vec2[T]) SetVec2(in *Vec3[T]) *Vec2[T] {
	v.X = in.X
	v.Y = in.Y
	return v
}

func (v *Vec2[T]) SetVec3(in *Vec3[T]) *Vec2[T] {
	v.X = in.X
	v.Y = in.Y

	return v
}

func (v *Vec2[T]) SetVec4(in *Vec4[T]) *Vec2[T] {
	v.X = in.X
	v.Y = in.Y

	return v
}

func (v *Vec2[T]) SetHomogenousVec3(in *Vec3[T]) *Vec2[T] {
	factor := T(1) / in.Z

	v.X = in.X * factor
	v.Y = in.Y * factor

	return v
}

func (v *Vec2[T]) SetHomogenousVec4(in *Vec4[T]) *Vec2[T] {
	factor := T(1) / in.W

	v.X = in.X * factor
	v.Y = in.Y * factor

	return v
}

func (v *Vec2[T]) Normalize() *Vec2[T] {
	vecLen := v.Len()

	if abs[T](vecLen-1) > 0.0001 {
		return v
	}

	inverse := 1.0 / vecLen
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

func (v *Vec2[T]) AddVec2(other *Vec2[T]) *Vec2[T] {
	v.X += other.X
	v.Y += other.Y

	return v
}

func (v *Vec2[T]) SubtractVec2(other *Vec2[T]) *Vec2[T] {
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

func (v *Vec2[T]) GreaterThan(other *Vec2[T]) bool {
	if v.X <= other.X || v.Y <= other.Y {
		return false
	}

	return true
}

func (v *Vec2[T]) GreaterThanEqual(other *Vec2[T]) bool {
	if v.X < other.X || v.Y < other.Y {
		return false
	}

	return true
}

func (v *Vec2[T]) LessThan(other *Vec2[T]) bool {
	if v.X >= other.X || v.Y >= other.Y {
		return false
	}

	return true
}

func (v *Vec2[T]) LessThanEqual(other *Vec2[T]) bool {
	if v.X > other.X || v.Y > other.Y {
		return false
	}

	return true
}

func (v *Vec2[T]) Lerp(other *Vec2[T], delta T) *Vec2[T] {
	v.X = v.X*(1-delta) + other.X*delta
	v.Y = v.Y*(1-delta) + other.Y*delta

	return v
}

func (v *Vec2[T]) SetLinearGradient(point0, point1, position *Vec2[T]) T {
	lineX := point1.X - point0.X
	lineY := point1.Y - point0.Y

	return (lineX*(position.X-point0.X) + lineY*(position.Y-point0.Y)) / (lineX*lineX + lineY*lineY)
}

func (v *Vec2[T]) Rotate(angleRad float64) *Vec2[T] {
	cos := T(math.Cos(angleRad))
	sin := T(math.Sin(angleRad))

	x := v.X*cos - v.Y*sin
	y := v.X*sin + v.Y*cos

	v.X = x
	v.Y = y
	return v
}

func (v *Vec2[T]) Transform(m *Mat2x2[T]) *Vec2[T] {
	x := m[0][0]*v.X + m[1][0]*v.Y
	y := m[0][1]*v.X + m[1][1]*v.Y

	v.X = x
	v.Y = y

	return v
}

func (v *Vec2[T]) TransformHomogenous(m *Mat3x3[T]) *Vec2[T] {
	x := m[0][0]*v.X + m[1][0]*v.Y + m[2][0]
	y := m[0][1]*v.X + m[1][1]*v.Y + m[2][1]
	z := m[0][2]*v.X + m[1][2]*v.Y + m[2][2]

	factor := 1.0 / z
	v.X = x * factor
	v.Y = y * factor

	return v
}
