package math

import "math"

type Vec2 struct {
	X float32
	Y float32
}

func (v *Vec2) SetVec2(in *Vec3) *Vec2 {
	v.X = in.X
	v.Y = in.Y
	return v
}

func (v *Vec2) SetVec3(in *Vec3) *Vec2 {
	v.X = in.X
	v.Y = in.Y

	return v
}

func (v *Vec2) SetVec4(in *Vec4) *Vec2 {
	v.X = in.X
	v.Y = in.Y

	return v
}

func (v *Vec2) SetHomogenousVec3(in *Vec3) *Vec2 {
	factor := 1 / in.Z

	v.X = in.X * factor
	v.Y = in.Y * factor

	return v
}

func (v *Vec2) SetHomogenousVec4(in *Vec4) *Vec2 {
	factor := 1 / in.W

	v.X = in.X * factor
	v.Y = in.Y * factor

	return v
}

func (v *Vec2) Normalize() *Vec2 {
	vecLen := v.Len()

	if abs(vecLen-1) > 0.0001 {
		return v
	}

	inverse := 1.0 / vecLen
	v.X *= inverse
	v.Y *= inverse

	return v
}

func (v *Vec2) Len() float32 {
	sqr := float64(v.X*v.X + v.Y*v.Y)
	return float32(math.Sqrt(sqr))
}

func (v *Vec2) LenSqr() float32 {
	return v.X*v.X + v.Y*v.Y
}

func (v *Vec2) AddVec2(other *Vec2) *Vec2 {
	v.X += other.X
	v.Y += other.Y

	return v
}

func (v *Vec2) SubtractVec2(other *Vec2) *Vec2 {
	v.X -= other.X
	v.Y -= other.Y

	return v
}

func (v *Vec2) DotProduct(other *Vec2) float32 {
	return v.X*other.X + v.Y*other.Y
}

func (v *Vec2) Scale(scale float32) *Vec2 {
	v.X *= scale
	v.Y *= scale

	return v
}

func (v *Vec2) Equal(other *Vec2, epsilon float32) bool {
	if abs(v.X-other.X) > epsilon {
		return false
	}

	if abs(v.Y-other.Y) > epsilon {
		return false
	}

	return true
}

func (v *Vec2) GreaterThan(other *Vec2) bool {
	if v.X <= other.X || v.Y <= other.Y {
		return false
	}

	return true
}

func (v *Vec2) GreaterThanEqual(other *Vec2) bool {
	if v.X < other.X || v.Y < other.Y {
		return false
	}

	return true
}

func (v *Vec2) LessThan(other *Vec2) bool {
	if v.X >= other.X || v.Y >= other.Y {
		return false
	}

	return true
}

func (v *Vec2) LessThanEqual(other *Vec2) bool {
	if v.X > other.X || v.Y > other.Y {
		return false
	}

	return true
}

func (v *Vec2) Lerp(other *Vec2, delta float32) *Vec2 {
	v.X = v.X*(1-delta) + other.X*delta
	v.Y = v.Y*(1-delta) + other.Y*delta

	return v
}

func (v *Vec2) SetLinearGradient(point0, point1, position *Vec2) float32 {
	lineX := point1.X - point0.X
	lineY := point1.Y - point0.Y

	return (lineX*(position.X-point0.X) + lineY*(position.Y-point0.Y)) / (lineX*lineX + lineY*lineY)
}

func (v *Vec2) Rotate(angleRad float32) *Vec2 {
	cos := float32(math.Cos(float64(angleRad)))
	sin := float32(math.Sin(float64(angleRad)))

	x := v.X*cos - v.Y*sin
	y := v.X*sin + v.Y*cos

	v.X = x
	v.Y = y
	return v
}

func (v *Vec2) Transform(m *Mat2x2) *Vec2 {
	x := m[0][0]*v.X + m[1][0]*v.Y
	y := m[0][1]*v.X + m[1][1]*v.Y

	v.X = x
	v.Y = y

	return v
}

func (v *Vec2) TransformHomogenous(m *Mat3x3) *Vec2 {
	x := m[0][0]*v.X + m[1][0]*v.Y + m[2][0]
	y := m[0][1]*v.X + m[1][1]*v.Y + m[2][1]
	z := m[0][2]*v.X + m[1][2]*v.Y + m[2][2]

	factor := 1.0 / z
	v.X = x * factor
	v.Y = y * factor

	return v
}
