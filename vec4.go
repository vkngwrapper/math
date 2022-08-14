package math

import "math"

type Vec4 struct {
	X float32
	Y float32
	Z float32
	W float32
}

func (v *Vec4) SetVec4(in *Vec4) *Vec4 {
	v.X = in.X
	v.Y = in.Y
	v.Z = in.Z
	v.W = in.W

	return v
}

func (v *Vec4) SetVec3(in *Vec3) *Vec4 {
	v.X = in.X
	v.Y = in.Y
	v.Z = in.Z
	v.W = 1

	return v
}

func (v *Vec4) SetVec2(in *Vec2) *Vec4 {
	v.X = in.X
	v.Y = in.Y
	v.Z = 0
	v.W = 1

	return v
}

func (v *Vec4) Normalize() *Vec4 {
	vecLen := v.Len()

	if math.Abs(float64(vecLen-1)) > 0.0001 {
		return v
	}

	inverse := 1.0 / vecLen
	v.X *= inverse
	v.Y *= inverse
	v.Z *= inverse
	v.W *= inverse

	return v
}

func (v *Vec4) Len() float32 {
	sqr := float64(v.X*v.X + v.Y*v.Y + v.Z*v.Z + v.W*v.W)
	return float32(math.Sqrt(sqr))
}

func (v *Vec4) LenSqr() float32 {
	return v.X*v.X + v.Y*v.Y + v.Z*v.Z + v.W*v.W
}

func (v *Vec4) AddVec4(other *Vec4) *Vec4 {
	v.X += other.X
	v.Y += other.Y
	v.Z += other.Z
	v.W += other.W

	return v
}

func (v *Vec4) SubtractVec4(other *Vec4) *Vec4 {
	v.X -= other.X
	v.Y -= other.Y
	v.Z -= other.Z
	v.W -= other.W

	return v
}

func (v *Vec4) DotProduct(other *Vec4) float32 {
	return v.X*other.X + v.Y*other.Y + v.Z*other.Z + v.W*other.W
}

func (v *Vec4) Scale(scale float32) *Vec4 {
	v.X *= scale
	v.Y *= scale
	v.Z *= scale
	v.W *= scale

	return v
}

func (v *Vec4) Equal(other *Vec4, epsilon float32) bool {
	if float32(math.Abs(float64(v.X-other.X))) > epsilon {
		return false
	}

	if abs(v.Y-other.Y) > epsilon {
		return false
	}

	if abs(v.Z-other.Z) > epsilon {
		return false
	}

	if abs(v.W-other.W) > epsilon {
		return false
	}

	return true
}

func (v *Vec4) GreaterThan(other *Vec4) bool {
	if v.X <= other.X || v.Y <= other.Y || v.Z <= other.Z || v.W <= other.W {
		return false
	}

	return true
}

func (v *Vec4) GreaterThanEqual(other *Vec4) bool {
	if v.X < other.X || v.Y < other.Y || v.Z < other.Z || v.W < other.W {
		return false
	}

	return true
}

func (v *Vec4) LessThan(other *Vec4) bool {
	if v.X >= other.X || v.Y >= other.Y || v.Z >= other.Z || v.W >= other.W {
		return false
	}

	return true
}

func (v *Vec4) LessThanEqual(other *Vec4) bool {
	if v.X > other.X || v.Y > other.Y || v.Z > other.Z || v.W > other.W {
		return false
	}

	return true
}

func (v *Vec4) Transform(m *Mat4x4) *Vec4 {
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

func (v *Vec4) Lerp(other *Vec4, delta float32) *Vec4 {
	v.X = v.X*(1-delta) + other.X*delta
	v.Y = v.Y*(1-delta) + other.Y*delta
	v.Z = v.Z*(1-delta) + other.Z*delta
	v.W = v.W*(1-delta) + other.W*delta

	return v
}

func (v *Vec4) Rotate(angleRad float32, normal *Vec3) *Vec4 {
	cos := float32(math.Cos(float64(angleRad)))
	sin := float32(math.Sin(float64(angleRad)))
	inverseCos := 1 - cos

	var unitAxis Vec3
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

func (v *Vec4) RotateWithQuaternion(q *Quaternion) *Vec4 {
	quatVector := Vec3{q.X, q.Y, q.Z}
	vec3 := Vec3{v.X, v.Y, v.Z}

	var uv Vec3
	uv.SetVec3(&quatVector).CrossProduct(&vec3)

	var uuv Vec3
	uuv.SetVec3(&quatVector).CrossProduct(&uv)

	v.X += ((uv.X * q.W) + uuv.X) * 2.0
	v.Y += ((uv.Y * q.W) + uuv.Y) * 2.0
	v.Z += ((uv.Z * q.W) + uuv.Z) * 2.0

	return v
}

func (v *Vec4) RotateX(angleRad float32) *Vec4 {
	cos := float32(math.Cos(float64(angleRad)))
	sin := float32(math.Sin(float64(angleRad)))

	y := v.Y*cos - v.Z*sin
	z := v.Y*sin + v.Z*cos

	v.Y = y
	v.Z = z

	return v
}

func (v *Vec4) RotateY(angleRad float32) *Vec4 {
	cos := float32(math.Cos(float64(angleRad)))
	sin := float32(math.Sin(float64(angleRad)))

	x := v.X*cos + v.Z*sin
	z := -v.X*sin + v.Z*cos

	v.X = x
	v.Z = z

	return v
}

func (v *Vec4) RotateZ(angleRad float32) *Vec4 {
	cos := float32(math.Cos(float64(angleRad)))
	sin := float32(math.Sin(float64(angleRad)))

	x := v.X*cos - v.Y*sin
	y := v.X*sin + v.Y*cos

	v.X = x
	v.Y = y

	return v
}
