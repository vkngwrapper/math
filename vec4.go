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

func (v *Vec4[T]) SetTransform(input *Vec4[T], m *Mat4x4[T]) *Vec4[T] {
	v.X = m[0][0]*input.X + m[1][0]*input.Y + m[2][0]*input.Z + m[3][0]*input.W
	v.Y = m[0][1]*input.X + m[1][1]*input.Y + m[2][1]*input.Z + m[3][1]*input.W
	v.Z = m[0][2]*input.X + m[1][2]*input.Y + m[2][2]*input.Z + m[3][2]*input.W
	v.W = m[0][3]*input.X + m[1][3]*input.Y + m[2][3]*input.Z + m[3][3]*input.W

	return v
}

func (v *Vec4[T]) SetTransformHomogenous(input *Vec4[T], m *Mat4x4[T]) *Vec4[T] {
	w := m[0][3]*input.X + m[1][3]*input.Y + m[2][3]*input.Z + m[3][3]*input.W

	v.X = (m[0][0]*input.X + m[1][0]*input.Y + m[2][0]*input.Z + m[3][0]*input.W) / w
	v.Y = (m[0][1]*input.X + m[1][1]*input.Y + m[2][1]*input.Z + m[3][1]*input.W) / w
	v.Z = (m[0][2]*input.X + m[1][2]*input.Y + m[2][2]*input.Z + m[3][2]*input.W) / w
	v.W = 1.0

	return v
}

func (v *Vec4[T]) SetRotateWithQuaternion(input *Vec4[T], q *Quaternion[T]) *Vec4[T] {
	quatVector := Vec3[T]{q.X, q.Y, q.Z}
	vec3 := Vec3[T]{input.X, input.Y, input.Z}
	var uv Vec3[T]

	uv.SetCrossProduct(&quatVector, &vec3)

	var uuv Vec3[T]
	uuv.SetCrossProduct(&quatVector, &uv)

	v.X = input.X + ((uv.X*q.W)+uuv.X)*2.0
	v.Y = input.Y + ((uv.Y*q.W)+uuv.Y)*2.0
	v.Z = input.Z + ((uv.Z*q.W)+uuv.Z)*2.0
	v.W = input.W

	return v
}

func (v *Vec4[T]) SetRotate(input *Vec4[T], angleRad float64, normal *Vec3[T]) *Vec4[T] {
	cos := T(math.Cos(angleRad))
	sin := T(math.Sin(angleRad))
	inverseCos := 1 - cos

	var unitAxis Vec3[T]
	unitAxis.SetVec3(normal).Normalize()

	v.X = (inverseCos*unitAxis.X*unitAxis.X+cos)*input.X +
		(inverseCos*unitAxis.X*unitAxis.Y+unitAxis.Z*sin)*input.Y +
		(inverseCos*unitAxis.X*unitAxis.Z-unitAxis.Y*sin)*input.Z
	v.Y = (inverseCos*unitAxis.X*unitAxis.Y-unitAxis.Z*sin)*input.X +
		(inverseCos*unitAxis.Y*unitAxis.Y+cos)*input.Y +
		(inverseCos*unitAxis.Y*unitAxis.Z+unitAxis.X*sin)*input.Z
	v.Z = (inverseCos*unitAxis.X*unitAxis.Z+unitAxis.Y*sin)*input.X +
		(inverseCos*unitAxis.Y*unitAxis.Z-unitAxis.X*sin)*input.Y +
		(inverseCos*unitAxis.Z*unitAxis.Z+cos)*input.Z
	v.W = input.W

	return v
}

func (v *Vec4[T]) SetRotateX(input *Vec4[T], angleRad float64) *Vec4[T] {
	cos := T(math.Cos(angleRad))
	sin := T(math.Sin(angleRad))

	y := input.Y*cos - input.Z*sin
	z := input.Y*sin + input.Z*cos

	v.X = input.X
	v.Y = y
	v.Z = z
	v.W = input.W

	return v
}

func (v *Vec4[T]) SetRotateY(input *Vec4[T], angleRad float64) *Vec4[T] {
	cos := T(math.Cos(angleRad))
	sin := T(math.Sin(angleRad))

	x := input.X*cos + input.Z*sin
	z := -input.X*sin + input.Z*cos

	v.X = x
	v.Y = input.Y
	v.Z = z
	v.W = input.W

	return v
}

func (v *Vec4[T]) SetRotateZ(input *Vec4[T], angleRad float64) *Vec4[T] {
	cos := T(math.Cos(angleRad))
	sin := T(math.Sin(angleRad))

	x := input.X*cos - input.Y*sin
	y := input.X*sin + input.Y*cos

	v.X = x
	v.Y = y
	v.Z = input.Z
	v.W = input.W

	return v
}

func (v *Vec4[T]) SetNormalizeVec4(input *Vec4[T]) *Vec4[T] {
	vecLen := input.Len()

	if abs[T](vecLen-1) < 0.0001 {
		v.X = input.X
		v.Y = input.Y
		v.Z = input.Z
		v.W = input.W
		return v
	}

	inverse := 1.0 / vecLen
	v.X = input.X * inverse
	v.Y = input.Y * inverse
	v.Z = input.Z * inverse
	v.W = input.W * inverse

	return v
}

func (v *Vec4[T]) SetAddVec4(lhs *Vec4[T], rhs *Vec4[T]) *Vec4[T] {
	v.X = lhs.X + rhs.X
	v.Y = lhs.Y + rhs.Y
	v.Z = lhs.Z + rhs.Z
	v.W = lhs.W + rhs.W

	return v
}

func (v *Vec4[T]) SetSubtractVec4(lhs *Vec4[T], rhs *Vec4[T]) *Vec4[T] {
	v.X = lhs.X - rhs.X
	v.Y = lhs.Y - rhs.Y
	v.Z = lhs.Z - rhs.Z
	v.W = lhs.W - rhs.W

	return v
}

func (v *Vec4[T]) SetScaleVec4(input *Vec4[T], scale T) *Vec4[T] {
	v.X = input.X * scale
	v.Y = input.Y * scale
	v.Z = input.Z * scale
	v.W = input.W * scale

	return v
}

func (v *Vec4[T]) SetLerpVec4(lhs *Vec4[T], rhs *Vec4[T], delta T) *Vec4[T] {
	v.X = lhs.X*(1-delta) + rhs.X*delta
	v.Y = lhs.Y*(1-delta) + rhs.Y*delta
	v.Z = lhs.Z*(1-delta) + rhs.Z*delta
	v.W = lhs.W*(1-delta) + rhs.W*delta

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

func (v *Vec4[T]) TransformHomogenous(m *Mat4x4[T]) *Vec4[T] {
	x := m[0][0]*v.X + m[1][0]*v.Y + m[2][0]*v.Z + m[3][0]*v.W
	y := m[0][1]*v.X + m[1][1]*v.Y + m[2][1]*v.Z + m[3][1]*v.W
	z := m[0][2]*v.X + m[1][2]*v.Y + m[2][2]*v.Z + m[3][2]*v.W
	w := m[0][3]*v.X + m[1][3]*v.Y + m[2][3]*v.Z + m[3][3]*v.W

	v.X = x / w
	v.Y = y / w
	v.Z = z / w
	v.W = 1.0

	return v
}

func (v *Vec4[T]) Lerp(other *Vec4[T], delta T) *Vec4[T] {
	v.X = v.X*(1-delta) + other.X*delta
	v.Y = v.Y*(1-delta) + other.Y*delta
	v.Z = v.Z*(1-delta) + other.Z*delta
	v.W = v.W*(1-delta) + other.W*delta

	return v
}

func (v *Vec4[T]) Rotate(angleRad float64, normal *Vec3[T]) *Vec4[T] {
	cos := T(math.Cos(angleRad))
	sin := T(math.Sin(angleRad))
	inverseCos := 1 - cos

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

func (v *Vec4[T]) RotateX(angleRad float64) *Vec4[T] {
	cos := T(math.Cos(angleRad))
	sin := T(math.Sin(angleRad))

	y := v.Y*cos - v.Z*sin
	z := v.Y*sin + v.Z*cos

	v.Y = y
	v.Z = z

	return v
}

func (v *Vec4[T]) RotateY(angleRad float64) *Vec4[T] {
	cos := T(math.Cos(angleRad))
	sin := T(math.Sin(angleRad))

	x := v.X*cos + v.Z*sin
	z := -v.X*sin + v.Z*cos

	v.X = x
	v.Z = z

	return v
}

func (v *Vec4[T]) RotateZ(angleRad float64) *Vec4[T] {
	cos := T(math.Cos(angleRad))
	sin := T(math.Sin(angleRad))

	x := v.X*cos - v.Y*sin
	y := v.X*sin + v.Y*cos

	v.X = x
	v.Y = y

	return v
}
