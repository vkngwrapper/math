package math

import "math"

type Vec4[T FloatingPoint] struct {
	X T
	Y T
	Z T
	W T
}

func (v *Vec4[T]) SetVec4(in *Vec4[T]) *Vec4[T] {
	v.X = in.X
	v.Y = in.Y
	v.Z = in.Z
	v.W = in.W

	return v
}

func (v *Vec4[T]) SetVec3(in *Vec3[T]) *Vec4[T] {
	v.X = in.X
	v.Y = in.Y
	v.Z = in.Z
	v.W = 1

	return v
}

func (v *Vec4[T]) SetVec2(in *Vec2[T]) *Vec4[T] {
	v.X = in.X
	v.Y = in.Y
	v.Z = 0
	v.W = 1

	return v
}

func (v *Vec4[T]) Normalize() *Vec4[T] {
	vecLen := v.Len()

	if abs[T](vecLen-1) > 0.0001 {
		return v
	}

	inverse := 1.0 / vecLen
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

func (v *Vec4[T]) GreaterThan(other *Vec4[T]) bool {
	if v.X <= other.X || v.Y <= other.Y || v.Z <= other.Z || v.W <= other.W {
		return false
	}

	return true
}

func (v *Vec4[T]) GreaterThanEqual(other *Vec4[T]) bool {
	if v.X < other.X || v.Y < other.Y || v.Z < other.Z || v.W < other.W {
		return false
	}

	return true
}

func (v *Vec4[T]) LessThan(other *Vec4[T]) bool {
	if v.X >= other.X || v.Y >= other.Y || v.Z >= other.Z || v.W >= other.W {
		return false
	}

	return true
}

func (v *Vec4[T]) LessThanEqual(other *Vec4[T]) bool {
	if v.X > other.X || v.Y > other.Y || v.Z > other.Z || v.W > other.W {
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

func (v *Vec4[T]) Lerp(other *Vec4[T], delta T) *Vec4[T] {
	v.X = v.X*(1-delta) + other.X*delta
	v.Y = v.Y*(1-delta) + other.Y*delta
	v.Z = v.Z*(1-delta) + other.Z*delta
	v.W = v.W*(1-delta) + other.W*delta

	return v
}

func (v *Vec4[T]) Rotate(angleRad T, normal *Vec3[T]) *Vec4[T] {
	cos := T(math.Cos(float64(angleRad)))
	sin := T(math.Sin(float64(angleRad)))
	inverseCos := T(1) - cos

	var unitAxis Vec3[T]
	unitAxis.SetVec3(normal).Normalize()

	x := (inverseCos*unitAxis.X*unitAxis.X+cos)*v.X +
		(inverseCos*unitAxis.X*unitAxis.Y+unitAxis.Z*sin)*v.Y +
		(inverseCos*unitAxis.X*unitAxis.Z-unitAxis.Y*sin)*v.Z

	y := (inverseCos*unitAxis.X*unitAxis.Y-unitAxis.Z*sin)*v.X +
		(inverseCos*unitAxis.Y*unitAxis.Y+cos)*v.Y +
		(inverseCos*unitAxis.Y*unitAxis.Z+unitAxis.X*sin)*v.Z

	z := (inverseCos*unitAxis.X*unitAxis.Z+unitAxis.Y*sin)*v.X +
		(inverseCos*unitAxis.Y*unitAxis.Z-unitAxis.X*sin)*v.Y +
		(inverseCos*unitAxis.Z*unitAxis.Z+cos)*v.Z

	v.X = x
	v.Y = y
	v.Z = z

	return v
}

func (v *Vec4[T]) RotateWithQuaternion(q *Quaternion[T]) *Vec4[T] {
	quatVector := Vec3[T]{q.X, q.Y, q.Z}
	vec3 := Vec3[T]{v.X, v.Y, v.Z}

	var uv Vec3[T]
	uv.SetCrossProduct(&quatVector, &vec3)

	var uuv Vec3[T]
	uuv.SetCrossProduct(&quatVector, &uv)

	v.X += ((uv.X * q.W) + uuv.X) * 2.0
	v.Y += ((uv.Y * q.W) + uuv.Y) * 2.0
	v.Z += ((uv.Z * q.W) + uuv.Z) * 2.0

	return v
}

func (v *Vec4[T]) RotateX(angleRad T) *Vec4[T] {
	cos := T(math.Cos(float64(angleRad)))
	sin := T(math.Sin(float64(angleRad)))

	y := v.Y*cos - v.Z*sin
	z := v.Y*sin + v.Z*cos

	v.Y = y
	v.Z = z

	return v
}

func (v *Vec4[T]) RotateY(angleRad T) *Vec4[T] {
	cos := T(math.Cos(float64(angleRad)))
	sin := T(math.Sin(float64(angleRad)))

	x := v.X*cos + v.Z*sin
	z := -v.X*sin + v.Z*cos

	v.X = x
	v.Z = z

	return v
}

func (v *Vec4[T]) RotateZ(angleRad T) *Vec4[T] {
	cos := T(math.Cos(float64(angleRad)))
	sin := T(math.Sin(float64(angleRad)))

	x := v.X*cos - v.Y*sin
	y := v.X*sin + v.Y*cos

	v.X = x
	v.Y = y

	return v
}
