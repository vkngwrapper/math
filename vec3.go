package math

import (
	"math"
)

type Vec3[T FloatingPoint] struct {
	X T
	Y T
	Z T
}

func (v *Vec3[T]) SetVec3(in *Vec3[T]) *Vec3[T] {
	v.X = in.X
	v.Y = in.Y
	v.Z = in.Z

	return v
}

func (v *Vec3[T]) SetVec2(in *Vec2[T]) *Vec3[T] {
	v.X = in.X
	v.Y = in.Y
	v.Z = 1

	return v
}

func (v *Vec3[T]) SetVec4(in *Vec4[T]) *Vec3[T] {
	v.X = in.X
	v.Y = in.Y
	v.Z = in.Z

	return v
}

func (v *Vec3[T]) SetHomogenousVec4(in *Vec4[T]) *Vec3[T] {
	factor := T(1) / in.W

	v.X = in.X * factor
	v.Y = in.Y * factor
	v.Z = in.Z * factor

	return v
}

func (v *Vec3[T]) SetCrossProduct(lhs, rhs *Vec3[T]) *Vec3[T] {
	v.X = lhs.Y*rhs.Z - lhs.Z*rhs.Y
	v.Y = lhs.Z*rhs.X - lhs.X*rhs.Z
	v.Z = lhs.X*rhs.Y - lhs.Y*rhs.X

	return v
}

func (v *Vec3[T]) Normalize() *Vec3[T] {
	vecLen := v.Len()

	if abs[T](vecLen-1) > 0.0001 {
		return v
	}

	inverse := 1.0 / vecLen
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

func (v *Vec3[T]) CrossProductQuaternion(other *Quaternion[T]) *Vec3[T] {
	quatVector := Vec3[T]{other.X, other.Y, other.Z}
	var uv, uuv Vec3[T]

	uv.SetVec3(&quatVector).CrossProduct(v)
	uuv.SetVec3(&quatVector).CrossProduct(&uv)

	v.X = v.X + ((uv.X*other.W)+uuv.X)*2.0
	v.Y = v.Y + ((uv.Y*other.W)+uuv.Y)*2.0
	v.Z = v.Z + ((uv.Z*other.W)+uuv.Z)*2.0

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

func (v *Vec3[T]) Equal(other *Vec3[T], epsilon T) bool {
	if abs[T](v.X-other.X) > epsilon {
		return false
	}

	if abs[T](v.Y-other.Y) > epsilon {
		return false
	}

	if abs[T](v.Z-other.Z) > epsilon {
		return false
	}

	return true
}

func (v *Vec3[T]) GreaterThan(other *Vec3[T]) bool {
	if v.X <= other.X || v.Y <= other.Y || v.Z <= other.Z {
		return false
	}

	return true
}

func (v *Vec3[T]) GreaterThanEqual(other *Vec3[T]) bool {
	if v.X < other.X || v.Y < other.Y || v.Z < other.Z {
		return false
	}

	return true
}

func (v *Vec3[T]) LessThan(other *Vec3[T]) bool {
	if v.X >= other.X || v.Y >= other.Y || v.Z >= other.Z {
		return false
	}

	return true
}

func (v *Vec3[T]) LessThanEqual(other *Vec3[T]) bool {
	if v.X > other.X || v.Y > other.Y || v.Z > other.Z {
		return false
	}

	return true
}

func (v *Vec3[T]) Lerp(other *Vec3[T], delta T) *Vec3[T] {
	v.X = v.X*(1-delta) + other.X*delta
	v.Y = v.Y*(1-delta) + other.Y*delta
	v.Z = v.Z*(1-delta) + other.Z*delta

	return v
}

func (v *Vec3[T]) SetClosestPointOnLine(point, endpoint1, endpoint2 *Vec3[T]) *Vec3[T] {

	var unitLine Vec3[T]
	unitLine.SetVec3(endpoint2).SubtractVec3(endpoint1)

	distanceSquared := unitLine.LenSqr()

	unitLine.Normalize()

	var lineToPoint Vec3[T]
	lineToPoint.SetVec3(point).SubtractVec3(endpoint1)

	distance := lineToPoint.DotProduct(&unitLine)

	if distance <= 0 {
		return endpoint1
	} else if distance*distance >= distanceSquared {
		return endpoint2
	}

	v.X = endpoint1.X + unitLine.X*distance
	v.Y = endpoint1.Y + unitLine.Y*distance
	v.Z = endpoint1.Z + unitLine.Z*distance

	return v
}

func (v *Vec3[T]) SetEuclidianFromPolarVec2(polar *Vec2[T]) *Vec3[T] {
	latitudeCos := T(math.Cos(float64(polar.X)))
	latitudeSin := T(math.Sin(float64(polar.X)))
	longitudeCos := T(math.Cos(float64(polar.Y)))
	longitudeSin := T(math.Sin(float64(polar.Y)))

	v.X = latitudeCos * longitudeSin
	v.Y = latitudeSin
	v.Z = latitudeCos * longitudeCos

	return v
}

func (v *Vec3[T]) SetPolarFromEuclidianVec3(euclidian *Vec3[T]) *Vec3[T] {
	var unitEuclidian Vec3[T]
	unitEuclidian.Normalize()

	xzDist := T(math.Sqrt(float64(unitEuclidian.X*unitEuclidian.X + unitEuclidian.Z*unitEuclidian.Z)))

	v.X = T(math.Asin(float64(unitEuclidian.Y)))
	v.Y = T(math.Atan2(float64(unitEuclidian.X), float64(unitEuclidian.Z)))
	v.Z = xzDist

	return v
}

func (v *Vec3[T]) Orthonormalize(other *Vec3[T]) *Vec3[T] {
	dotProduct := v.DotProduct(other)

	v.X = v.X - other.X*dotProduct
	v.Y = v.Y - other.Y*dotProduct
	v.Z = v.Z - other.Z*dotProduct

	return v.Normalize()
}

func (v *Vec3[T]) Project(model *Mat4x4[T], proj *Mat4x4[T], viewport *Vec4[T]) *Vec3[T] {
	xModel := model[0][0]*v.X + model[0][1]*v.Y + model[0][2]*v.Z + model[0][3]
	yModel := model[1][0]*v.X + model[1][1]*v.Y + model[1][2]*v.Z + model[1][3]
	zModel := model[2][0]*v.X + model[2][1]*v.Y + model[2][2]*v.Z + model[2][3]
	wModel := model[3][0] + v.X + model[3][1]*v.Y + model[3][2]*v.Z + model[3][3]

	x := proj[0][0]*xModel + proj[0][1]*yModel + proj[0][2]*zModel + proj[0][3]*wModel
	y := proj[1][0]*xModel + proj[1][1]*yModel + proj[1][2]*zModel + proj[1][3]*wModel
	z := proj[2][0]*xModel + proj[2][1]*yModel + proj[2][2]*zModel + proj[2][3]*wModel
	w := proj[3][0]*xModel + proj[3][1]*yModel + proj[3][2]*zModel + proj[3][3]*wModel

	factor := 1 / w
	v.X = (x*factor*0.5+0.5)*viewport.Z + viewport.X
	v.Y = (y*factor*0.5+0.5)*viewport.W + viewport.Y
	v.Z = z * factor

	return v
}

func (v *Vec3[T]) Unproject(model *Mat4x4[T], proj *Mat4x4[T], viewport *Vec4[T]) *Vec3[T] {
	var inverse Mat4x4[T]
	inverse.SetMat4x4(proj).MultMatrix4x4(model).Inverse()

	x := (v.X - viewport.X) / viewport.Z
	y := (v.Y - viewport.Y) / viewport.W
	x = x*2.0 - 1.0
	y = y*2.0 - 1.0
	z := v.Z

	v.X = inverse[0][0]*x + inverse[0][1]*y + inverse[0][2]*z + inverse[0][3]
	v.Y = inverse[1][0]*x + inverse[1][1]*y + inverse[1][2]*z + inverse[1][3]
	v.Z = inverse[2][0]*x + inverse[2][1]*y + inverse[2][2]*z + inverse[2][3]
	w := inverse[3][0]*x + inverse[3][1]*y + inverse[3][2]*z + inverse[3][3]

	factor := 1.0 / w
	v.X *= factor
	v.Y *= factor
	v.Z *= factor

	return v
}

func (v *Vec3[T]) Rotate(angleRad T, normal *Vec3[T]) *Vec3[T] {
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

func (v *Vec3[T]) RotateWithQuaternion(q *Quaternion[T]) *Vec3[T] {
	quatVector := Vec3[T]{q.X, q.Y, q.Z}

	var uv Vec3[T]
	uv.SetVec3(&quatVector).CrossProduct(v)

	var uuv Vec3[T]
	uuv.SetVec3(&quatVector).CrossProduct(&uv)

	v.X += ((uv.X * q.W) + uuv.X) * 2.0
	v.Y += ((uv.Y * q.W) + uuv.Y) * 2.0
	v.Z += ((uv.Z * q.W) + uuv.Z) * 2.0

	return v
}

func (v *Vec3[T]) RotateX(angleRad T) *Vec3[T] {
	cos := T(math.Cos(float64(angleRad)))
	sin := T(math.Sin(float64(angleRad)))

	y := v.Y*cos - v.Z*sin
	z := v.Y*sin + v.Z*cos

	v.Y = y
	v.Z = z

	return v
}

func (v *Vec3[T]) RotateY(angleRad T) *Vec3[T] {
	cos := T(math.Cos(float64(angleRad)))
	sin := T(math.Sin(float64(angleRad)))

	x := v.X*cos + v.Z*sin
	z := -v.X*sin + v.Z*cos

	v.X = x
	v.Z = z

	return v
}

func (v *Vec3[T]) RotateZ(angleRad T) *Vec3[T] {
	cos := T(math.Cos(float64(angleRad)))
	sin := T(math.Sin(float64(angleRad)))

	x := v.X*cos - v.Y*sin
	y := v.X*sin + v.Y*cos

	v.X = x
	v.Y = y

	return v
}

func (v *Vec3[T]) SetTriangleNormal(point1, point2, point3 *Vec3[T]) *Vec3[T] {
	var edge1, edge2 Vec3[T]

	edge2.SetVec3(point1).SubtractVec3(point3)
	normal := edge1.SetVec3(point1).SubtractVec3(point2).CrossProduct(&edge2).Normalize()

	v.X = normal.X
	v.Y = normal.Y
	v.Z = normal.Z

	return v
}

func (v *Vec3[T]) Transform(m *Mat3x3[T]) *Vec3[T] {
	x := m[0][0]*v.X + m[1][0]*v.Y + m[2][0]*v.Z
	y := m[0][1]*v.X + m[1][1]*v.Y + m[2][1]*v.Z
	z := m[0][2]*v.X + m[1][2]*v.Y + m[2][2]*v.Z

	v.X = x
	v.Y = y
	v.Z = z

	return v
}

func (v *Vec3[T]) TransformHomogenous(m *Mat4x4[T]) *Vec3[T] {
	x := m[0][0]*v.X + m[1][0]*v.Y + m[2][0]*v.Z + m[3][0]
	y := m[0][1]*v.X + m[1][1]*v.Y + m[2][1]*v.Z + m[3][1]
	z := m[0][2]*v.X + m[1][2]*v.Y + m[2][2]*v.Z + m[3][2]
	w := m[0][3]*v.X + m[1][3]*v.Y + m[2][3]*v.Z + m[3][3]

	factor := 1.0 / w
	v.X = x * factor
	v.Y = y * factor
	v.Z = z * factor

	return v
}
