package math

import "math"

type Mat4Row[T FloatingPoint] [4]T
type Mat4x4[T FloatingPoint] [4]Mat4Row[T]

func (m *Mat4x4[T]) SetIdentity() *Mat4x4[T] {
	m[0][0] = 1
	m[0][1] = 0
	m[0][2] = 0
	m[0][3] = 0
	m[1][0] = 0
	m[1][1] = 1
	m[1][2] = 0
	m[1][3] = 0
	m[2][0] = 0
	m[2][1] = 0
	m[2][2] = 1
	m[2][3] = 0
	m[3][0] = 0
	m[3][1] = 0
	m[3][2] = 0
	m[3][3] = 1

	return m
}

func (m *Mat4x4[T]) SetDiagonalScalar(s T) *Mat4x4[T] {
	m[0][0] = s
	m[0][1] = 0
	m[0][2] = 0
	m[0][3] = 0
	m[1][0] = 0
	m[1][1] = s
	m[1][2] = 0
	m[1][3] = 0
	m[2][0] = 0
	m[2][1] = 0
	m[2][2] = s
	m[2][3] = 0
	m[3][0] = 0
	m[3][1] = 0
	m[3][2] = 0
	m[3][3] = s

	return m
}

func (m *Mat4x4[T]) SetDiagonalVector(v *Vec4[T]) *Mat4x4[T] {
	m[0][0] = v.X
	m[0][1] = 0
	m[0][2] = 0
	m[0][3] = 0
	m[1][0] = 0
	m[1][1] = v.Y
	m[1][2] = 0
	m[1][3] = 0
	m[2][0] = 0
	m[2][1] = 0
	m[2][2] = v.Z
	m[2][3] = 0
	m[3][0] = 0
	m[3][1] = 0
	m[3][2] = 0
	m[3][3] = v.W

	return m
}

func (m *Mat4x4[T]) SetFrustum(left, right, bottom, top, nearVal, farVal T) *Mat4x4[T] {
	m[0][0] = (T(2.0) * nearVal) / (right - left)
	m[0][1] = 0
	m[0][2] = 0
	m[0][3] = 0
	m[1][0] = 0
	m[1][1] = (T(2) * nearVal) / (top - bottom)
	m[1][2] = 0
	m[1][3] = 0
	m[2][0] = (right + left) / (right - left)
	m[2][1] = (top + bottom) / (top - bottom)
	m[2][2] = farVal / (nearVal - farVal)
	m[2][3] = T(-1)
	m[3][0] = 0
	m[3][1] = 0
	m[3][2] = -(farVal * nearVal) / (farVal - nearVal)
	m[3][3] = 0

	return m
}

func (m *Mat4x4[T]) SetAxisAngle(axis *Vec3[T], angleRad T) *Mat4x4[T] {
	cos := T(math.Cos(float64(angleRad)))
	sin := T(math.Sin(float64(angleRad)))
	inverseCos := T(1) - cos

	var unitAxis Vec3[T]
	unitAxis.SetVec3(axis).Normalize()

	m[0][0] = inverseCos*unitAxis.X*unitAxis.X + cos
	m[0][1] = inverseCos*unitAxis.X*unitAxis.Y + unitAxis.Z*sin
	m[0][2] = inverseCos*unitAxis.X*unitAxis.Z - unitAxis.Y*sin
	m[0][3] = 0
	m[1][0] = inverseCos*unitAxis.X*unitAxis.Y - unitAxis.Z*sin
	m[1][1] = inverseCos*unitAxis.Y*unitAxis.Y + cos
	m[1][2] = inverseCos*unitAxis.Y*unitAxis.Z + unitAxis.X*sin
	m[1][3] = 0
	m[2][0] = inverseCos*unitAxis.X*unitAxis.Z + unitAxis.Y*sin
	m[2][1] = inverseCos*unitAxis.Y*unitAxis.Z - unitAxis.X*sin
	m[2][2] = inverseCos*unitAxis.Z*unitAxis.Z + cos
	m[2][3] = 0
	m[3][0] = 0
	m[3][1] = 0
	m[3][2] = 0
	m[3][3] = 1

	return m
}

func (m *Mat4x4[T]) SetOrientation(normal *Vec3[T], up *Vec3[T]) *Mat4x4[T] {
	if normal.Equal(up, 0.0001) {
		m.SetIdentity()
		return m
	}

	var rotationAxis Vec3[T]
	rotationAxis.SetVec3(up).CrossProduct(normal)

	angle := T(math.Acos(float64(normal.DotProduct(up))))

	return m.SetAxisAngle(&rotationAxis, angle)
}

func (m *Mat4x4[T]) SetTranslation(x, y, z T) *Mat4x4[T] {
	m[0][0] = 1
	m[0][1] = 0
	m[0][2] = 0
	m[0][3] = 0
	m[1][0] = 0
	m[1][1] = 1
	m[1][2] = 0
	m[1][3] = 0
	m[2][0] = 0
	m[2][1] = 0
	m[2][2] = 1
	m[2][3] = 0
	m[3][0] = x
	m[3][1] = y
	m[3][2] = z
	m[3][3] = 1

	return m
}

func (m *Mat4x4[T]) SetScale(x, y, z T) *Mat4x4[T] {
	m[0][0] = x
	m[0][1] = 0
	m[0][2] = 0
	m[0][3] = 0
	m[1][0] = 0
	m[1][1] = y
	m[1][2] = 0
	m[1][3] = 0
	m[2][0] = 0
	m[2][1] = 0
	m[2][2] = z
	m[2][3] = 0
	m[3][0] = 0
	m[3][1] = 0
	m[3][2] = 0
	m[3][3] = 1

	return m
}

func (m *Mat4x4[T]) SetMat4x4(other *Mat4x4[T]) *Mat4x4[T] {
	m[0][0] = other[0][0]
	m[0][1] = other[0][1]
	m[0][2] = other[0][2]
	m[0][3] = other[0][3]
	m[1][0] = other[1][0]
	m[1][1] = other[1][1]
	m[1][2] = other[1][2]
	m[1][3] = other[1][3]
	m[2][0] = other[2][0]
	m[2][1] = other[2][1]
	m[2][2] = other[2][2]
	m[2][3] = other[2][3]
	m[3][0] = other[3][0]
	m[3][1] = other[3][1]
	m[3][2] = other[3][2]
	m[3][3] = other[3][3]

	return m
}

func (m *Mat4x4[T]) SetMat3x3(other *Mat3x3[T]) *Mat4x4[T] {
	m[0][0] = other[0][0]
	m[0][1] = other[0][1]
	m[0][2] = other[0][2]
	m[0][3] = 0
	m[1][0] = other[1][0]
	m[1][1] = other[1][1]
	m[1][2] = other[1][2]
	m[1][3] = 0
	m[2][0] = other[2][0]
	m[2][1] = other[2][1]
	m[2][2] = other[2][2]
	m[2][3] = 0
	m[3][0] = 0
	m[3][1] = 0
	m[3][2] = 0
	m[3][3] = 1

	return m
}

func (m *Mat4x4[T]) SetColMajor(c1, c2, c3, c4 *Vec4[T]) *Mat4x4[T] {
	m[0][0] = c1.X
	m[0][1] = c1.Y
	m[0][2] = c1.Z
	m[0][3] = c1.W
	m[1][0] = c2.X
	m[1][1] = c2.Y
	m[1][2] = c2.Z
	m[1][3] = c2.W
	m[2][0] = c3.X
	m[2][1] = c3.Y
	m[2][2] = c3.Z
	m[2][3] = c3.W
	m[3][0] = c4.X
	m[3][1] = c4.Y
	m[3][2] = c4.Z
	m[3][3] = c4.W

	return m
}

func (m *Mat4x4[T]) SetRowMajor(r1, r2, r3, r4 *Vec4[T]) *Mat4x4[T] {
	m[0][0] = r1.X
	m[1][0] = r1.Y
	m[2][0] = r2.Z
	m[3][0] = r3.W
	m[0][1] = r2.X
	m[1][1] = r2.Y
	m[2][1] = r2.Z
	m[3][1] = r2.W
	m[0][2] = r3.X
	m[1][2] = r3.Y
	m[2][2] = r3.Z
	m[3][2] = r3.W
	m[0][3] = r4.X
	m[1][3] = r4.Y
	m[2][3] = r4.Z
	m[3][3] = r4.W

	return m
}

func (m *Mat4x4[T]) SetEulerAngles(yawRad, pitchRad, rollRad T) *Mat4x4[T] {
	yawCos := T(math.Cos(float64(yawRad)))
	pitchCos := T(math.Cos(float64(pitchRad)))
	rollCos := T(math.Cos(float64(rollRad)))
	yawSin := T(math.Sin(float64(yawRad)))
	pitchSin := T(math.Sin(float64(pitchRad)))
	rollSin := T(math.Sin(float64(rollRad)))

	m[0][0] = yawCos*rollCos + yawSin*pitchSin*rollSin
	m[0][1] = rollSin * pitchSin
	m[0][2] = -yawSin*rollCos + yawCos*pitchSin*rollSin
	m[0][3] = 0
	m[1][0] = -yawCos*rollSin + yawSin*pitchSin*rollCos
	m[1][1] = rollCos * pitchCos
	m[1][2] = rollSin*yawSin + yawCos*pitchSin*rollCos
	m[1][3] = 0
	m[2][0] = yawSin * pitchCos
	m[2][1] = -pitchSin
	m[2][2] = yawCos * pitchCos
	m[2][3] = 0
	m[3][0] = 0
	m[3][1] = 0
	m[3][2] = 0
	m[3][3] = 1

	return m
}

func (m *Mat4x4[T]) SetMatrixRotationFrom(other *Mat4x4[T]) *Mat4x4[T] {
	m[0][0] = other[0][0]
	m[0][1] = other[0][1]
	m[0][2] = other[0][2]
	m[0][3] = T(0)
	m[1][0] = other[1][0]
	m[1][1] = other[1][1]
	m[1][2] = other[1][2]
	m[1][3] = T(0)
	m[2][0] = other[2][0]
	m[2][1] = other[2][1]
	m[2][2] = other[2][2]
	m[2][3] = T(0)
	m[3][0] = T(0)
	m[3][1] = T(0)
	m[3][2] = T(0)
	m[3][3] = T(1)

	return m
}

func (m *Mat4x4[T]) SetInfinitePerspective(fovY, aspectRatio, zNear T) *Mat4x4[T] {
	r := T(math.Tan(float64(fovY/T(2)))) * zNear
	left := -r * aspectRatio
	right := r * aspectRatio
	bottom := -r
	top := r

	m[0][0] = (T(2) * zNear) / (right - left)
	m[0][1] = 0
	m[0][2] = 0
	m[0][3] = 0
	m[1][0] = 0
	m[1][1] = (T(2) * zNear) / (top - bottom)
	m[1][2] = 0
	m[1][3] = 0
	m[2][0] = 0
	m[2][1] = 0
	m[2][2] = -T(1)
	m[2][3] = -T(1)
	m[3][0] = 0
	m[3][1] = 0
	m[3][2] = -T(2) * zNear
	m[3][3] = 0

	return m
}

func (m *Mat4x4[T]) SetOrthographic2D(left, right, bottom, top T) *Mat4x4[T] {
	m[0][0] = T(2) / (right - left)
	m[0][1] = 0
	m[0][2] = 0
	m[0][3] = 0
	m[1][0] = 0
	m[1][1] = T(2) / (top - bottom)
	m[1][2] = 0
	m[1][3] = 0
	m[2][0] = 0
	m[2][1] = 0
	m[2][2] = T(1)
	m[2][3] = 0
	m[3][0] = -(right + left) / (right - left)
	m[3][1] = -(top + bottom) / (top - bottom)
	m[3][2] = 0
	m[3][3] = 1

	return m
}

func (m *Mat4x4[T]) SetOrthographic(left, right, bottom, top, zNear, zFar T) *Mat4x4[T] {
	m[0][0] = T(2) / (right - left)
	m[0][1] = 0
	m[0][2] = 0
	m[0][3] = 0
	m[1][0] = 0
	m[1][1] = T(2) / (top - bottom)
	m[1][2] = 0
	m[1][3] = 0
	m[2][0] = 0
	m[2][1] = 0
	m[2][2] = -T(1) / (zFar - zNear)
	m[2][3] = 0
	m[3][0] = -(right + left) / (right - left)
	m[3][1] = -(top + bottom) / (top - bottom)
	m[3][2] = -zNear / (zFar - zNear)
	m[3][3] = 1

	return m
}

func (m *Mat4x4[T]) SetPerspective(fovY, aspectRatio, zNear, zFar T) *Mat4x4[T] {
	tanHalfFovY := T(math.Tan(float64(fovY / T(2))))

	m[0][0] = T(1) / (aspectRatio * tanHalfFovY)
	m[0][1] = 0
	m[0][2] = 0
	m[0][3] = 0
	m[1][0] = 0
	m[1][1] = T(1) / tanHalfFovY
	m[1][2] = 0
	m[1][3] = 0
	m[2][0] = 0
	m[2][1] = 0
	m[2][2] = zFar / (zNear - zFar)
	m[2][3] = T(-1)
	m[3][0] = 0
	m[3][1] = 0
	m[3][2] = -(zFar * zNear) / (zFar - zNear)
	m[3][3] = 0

	return m
}

func (m *Mat4x4[T]) SetPerspectiveFOV(fov, width, height, zNear, zFar T) *Mat4x4[T] {
	h := T(math.Cos(0.5*float64(fov))) / T(math.Sin(0.5*float64(fov)))
	w := h * height / width

	m[0][0] = w
	m[0][1] = 0
	m[0][2] = 0
	m[0][3] = 0
	m[1][0] = 0
	m[1][1] = h
	m[1][2] = 0
	m[1][3] = 0
	m[2][0] = 0
	m[2][1] = 0
	m[2][2] = zFar / (zNear - zFar)
	m[2][3] = T(-1)
	m[3][0] = 0
	m[3][1] = 0
	m[3][2] = -(zFar * zNear) / (zFar - zNear)
	m[3][3] = 0

	return m
}

func (m *Mat4x4[T]) SetPickMatrix(center *Vec2[T], delta *Vec2[T], viewport Vec4[T]) *Mat4x4[T] {
	m.SetIdentity()

	if !(delta.X > 0 && delta.Y > 0) {
		return m
	}

	var temp Vec3[T]
	temp.X = (viewport.Z - 2*(center.X-viewport.X)) / delta.X
	temp.Y = (viewport.W - 2*(center.Y-viewport.Y)) / delta.Y
	temp.Z = 0

	var scale Vec3[T]
	scale.X = viewport.Z / delta.X
	scale.Y = viewport.W / delta.Y
	scale.Z = 1

	return m.Translate(&temp).Scale(&scale)
}

func (m *Mat4x4[T]) SetLookAt(eyePosition *Vec3[T], target *Vec3[T], up *Vec3[T]) *Mat4x4[T] {
	var f Vec3[T]
	f.SetVec3(target).SubtractVec3(eyePosition).Normalize()

	var s Vec3[T]
	s.SetVec3(&f).CrossProduct(up).Normalize()

	var u Vec3[T]
	u.SetVec3(&s).CrossProduct(&f)

	m[0][0] = s.X
	m[1][0] = s.Y
	m[2][0] = s.Z
	m[3][0] = 0
	m[0][1] = u.X
	m[1][1] = u.Y
	m[2][1] = u.Z
	m[3][1] = 0
	m[0][2] = -f.X
	m[1][2] = -f.Y
	m[2][2] = -f.Z
	m[3][2] = 0
	m[3][0] = -s.DotProduct(eyePosition)
	m[3][1] = -u.DotProduct(eyePosition)
	m[3][2] = f.DotProduct(eyePosition)
	m[3][3] = 1

	return m
}

func (m *Mat4x4[T]) SetMatrixCrossProduct(vec *Vec3[T]) *Mat4x4[T] {
	m[0][0] = 0
	m[0][1] = vec.Z
	m[0][2] = -vec.Y
	m[0][3] = 0
	m[1][0] = -vec.Z
	m[1][1] = 0
	m[1][2] = vec.X
	m[1][3] = 0
	m[2][0] = vec.Y
	m[2][1] = -vec.X
	m[2][2] = 0
	m[2][3] = 0
	m[3][0] = 0
	m[3][1] = 0
	m[3][2] = 0
	m[3][3] = 0

	return m
}

func (m *Mat4x4[T]) SetMultMatrix4x4(lhs, rhs *Mat4x4[T]) *Mat4x4[T] {
	m00 := lhs[0][0]*rhs[0][0] + lhs[1][0]*rhs[0][1] + lhs[2][0]*rhs[0][2] + lhs[3][0]*rhs[0][3]
	m10 := lhs[0][0]*rhs[1][0] + lhs[1][0]*rhs[1][1] + lhs[2][0]*rhs[1][2] + lhs[3][0]*rhs[1][3]
	m20 := lhs[0][0]*rhs[2][0] + lhs[1][0]*rhs[2][1] + lhs[2][0]*rhs[2][2] + lhs[3][0]*rhs[2][3]
	m30 := lhs[0][0]*rhs[3][0] + lhs[1][0]*rhs[3][1] + lhs[2][0]*rhs[3][2] + lhs[3][0]*rhs[3][3]

	m01 := lhs[0][1]*rhs[0][0] + lhs[1][1]*rhs[0][1] + lhs[2][1]*rhs[0][2] + lhs[3][1]*rhs[0][3]
	m11 := lhs[0][1]*rhs[1][0] + lhs[1][1]*rhs[1][1] + lhs[2][1]*rhs[1][2] + lhs[3][1]*rhs[1][3]
	m21 := lhs[0][1]*rhs[2][0] + lhs[1][1]*rhs[2][1] + lhs[2][1]*rhs[2][2] + lhs[3][1]*rhs[2][3]
	m31 := lhs[0][1]*rhs[3][0] + lhs[1][1]*rhs[3][1] + lhs[2][1]*rhs[3][2] + lhs[3][1]*rhs[3][3]

	m02 := lhs[0][2]*rhs[0][0] + lhs[1][2]*rhs[0][1] + lhs[2][2]*rhs[0][2] + lhs[3][2]*rhs[0][3]
	m12 := lhs[0][2]*rhs[1][0] + lhs[1][2]*rhs[1][1] + lhs[2][2]*rhs[1][2] + lhs[3][2]*rhs[1][3]
	m22 := lhs[0][2]*rhs[2][0] + lhs[1][2]*rhs[2][1] + lhs[2][2]*rhs[2][2] + lhs[3][2]*rhs[2][3]
	m32 := lhs[0][2]*rhs[3][0] + lhs[1][2]*rhs[3][1] + lhs[2][2]*rhs[3][2] + lhs[3][2]*rhs[3][3]

	m03 := lhs[0][3]*rhs[0][0] + lhs[1][3]*rhs[0][1] + lhs[2][3]*rhs[0][2] + lhs[3][3]*rhs[0][3]
	m13 := lhs[0][3]*rhs[1][0] + lhs[1][3]*rhs[1][1] + lhs[2][3]*rhs[1][2] + lhs[3][3]*rhs[1][3]
	m23 := lhs[0][3]*rhs[2][0] + lhs[1][3]*rhs[2][1] + lhs[2][3]*rhs[2][2] + lhs[3][3]*rhs[2][3]
	m33 := lhs[0][3]*rhs[3][0] + lhs[1][3]*rhs[3][1] + lhs[2][3]*rhs[3][2] + lhs[3][3]*rhs[3][3]

	m[0][0] = m00
	m[0][1] = m01
	m[0][2] = m02
	m[0][3] = m03
	m[1][0] = m10
	m[1][1] = m11
	m[1][2] = m12
	m[1][3] = m13
	m[2][0] = m20
	m[2][1] = m21
	m[2][2] = m22
	m[2][3] = m23
	m[3][0] = m30
	m[3][1] = m31
	m[3][2] = m32
	m[3][3] = m33

	return m
}

func (m *Mat4x4[T]) SetInterpolate(lhs, rhs *Mat4x4[T], delta T) *Mat4x4[T] {
	var oldMatrix Mat4x4[T]
	oldMatrix.SetMat4x4(lhs)

	var thisRotation Mat4x4[T]
	thisRotation.SetMatrixRotationFrom(lhs)

	var transposedRotation Mat4x4[T]
	transposedRotation.SetMat4x4(&thisRotation).Transpose()

	var deltaRotation Mat4x4[T]
	deltaRotation.SetMultMatrix4x4(rhs, &transposedRotation)

	var deltaAxis Vec3[T]
	var deltaAngle T
	deltaRotation.GetAxisAngle(&deltaAxis, &deltaAngle)

	m.SetAxisAngle(&deltaAxis, deltaAngle).MultMatrix4x4(&thisRotation)
	m[3][0] = oldMatrix[3][0] + delta*(rhs[3][0]-oldMatrix[3][0])
	m[3][1] = oldMatrix[3][1] + delta*(rhs[3][1]-oldMatrix[3][1])
	m[3][2] = oldMatrix[3][2] + delta*(rhs[3][2]-oldMatrix[3][2])

	return m
}

func (m *Mat4x4[T]) GetAxisAngle(outAxis *Vec3[T], outAngleRad *T) {
	epsilon := T(0.01)
	epsilon2 := T(0.1)

	if (abs[T](m[1][0]-m[0][1]) < epsilon) &&
		(abs[T](m[2][0]-m[0][2]) < epsilon) &&
		(abs[T](m[2][1]-m[1][2]) < epsilon) {

		if (abs[T](m[1][0]+m[0][1]) < epsilon2) &&
			(abs[T](m[2][0]+m[0][2]) < epsilon2) &&
			(abs[T](m[2][1]+m[1][2]) < epsilon2) &&
			(abs[T](m[0][0]+m[1][1]+m[2][2]-T(3)) < epsilon2) {

			*outAngleRad = T(0)
			outAxis.X = T(1)
			outAxis.Y = T(0)
			outAxis.Z = T(0)
			return
		}

		*outAngleRad = T(math.Pi)
		xx := (m[0][0] + T(1)) * T(0.5)
		yy := (m[1][1] + T(1)) * T(0.5)
		zz := (m[2][2] + T(1)) * T(0.5)
		xy := (m[1][0] + m[0][1]) * T(0.25)
		xz := (m[2][0] + m[0][2]) * T(0.25)
		yz := (m[2][1] + m[1][2]) * T(0.25)

		if xx > yy && xx > zz {
			if xx < epsilon {
				outAxis.X = T(0)
				outAxis.Y = T(0.7071)
				outAxis.Z = T(0.7071)
			} else {
				outAxis.X = T(math.Sqrt(float64(xx)))
				outAxis.Y = xy / outAxis.X
				outAxis.Z = xz / outAxis.X
			}
		} else if yy > zz {
			if yy < epsilon {
				outAxis.X = T(0.7071)
				outAxis.Y = T(0)
				outAxis.Z = T(0.7071)
			} else {
				outAxis.Y = T(math.Sqrt(float64(yy)))
				outAxis.X = xy / outAxis.Y
				outAxis.Z = yz / outAxis.Z
			}
		} else {
			if zz < epsilon {
				outAxis.X = T(0.7071)
				outAxis.Y = T(0.7071)
				outAxis.Z = T(0)
			} else {
				outAxis.Z = T(math.Sqrt(float64(zz)))
				outAxis.X = xz / outAxis.Z
				outAxis.Y = yz / outAxis.Z
			}
		}
		return
	}

	sSquared := (m[2][1]-m[1][2])*(m[2][1]-m[1][2]) +
		(m[2][0]-m[0][2])*(m[2][0]-m[0][2]) +
		(m[1][0]-m[0][1])*(m[1][0]-m[0][1])
	s := T(math.Sqrt(float64(sSquared)))

	if s < T(0.001) {
		s = T(1)
	}

	angleCos := (m[0][0] + m[1][1] + m[2][2] - T(1)) * T(0.5)
	if angleCos-T(1) < epsilon {
		*outAngleRad = T(math.Pi) * 0.25
	} else {
		*outAngleRad = T(math.Acos(float64(angleCos)))
	}

	outAxis.X = (m[1][2] - m[2][1]) / s
	outAxis.Y = (m[2][0] - m[0][2]) / s
	outAxis.Z = (m[0][1] - m[1][0]) / s
}

func (m *Mat4x4[T]) Transpose() *Mat4x4[T] {
	m[1][0], m[0][1] = m[0][1], m[1][0]
	m[2][0], m[0][2] = m[0][2], m[2][0]
	m[2][1], m[1][2] = m[1][2], m[2][1]

	m[3][0], m[0][3] = m[0][3], m[3][0]
	m[3][1], m[1][3] = m[1][3], m[3][1]
	m[3][2], m[2][3] = m[2][3], m[3][2]

	return m
}

func (m *Mat4x4[T]) Inverse() *Mat4x4[T] {
	coefficient0 := m[2][2]*m[3][3] - m[3][2]*m[2][3]
	coefficient2 := m[1][2]*m[3][3] - m[3][2]*m[1][3]
	coefficient3 := m[1][2]*m[2][3] - m[2][2]*m[1][3]

	coefficient4 := m[2][1]*m[3][3] - m[3][1]*m[2][3]
	coefficient6 := m[1][1]*m[3][3] - m[3][1]*m[1][3]
	coefficient7 := m[1][1]*m[2][3] - m[2][1]*m[1][3]

	coefficient8 := m[2][1]*m[3][2] - m[3][1]*m[2][2]
	coefficient10 := m[1][1]*m[3][2] - m[3][1]*m[1][2]
	coefficient11 := m[1][1]*m[2][2] - m[2][1]*m[1][2]

	coefficient12 := m[2][0]*m[3][3] - m[3][0]*m[2][3]
	coefficient14 := m[1][0]*m[3][3] - m[3][0]*m[1][3]
	coefficient15 := m[1][0]*m[2][3] - m[2][0]*m[1][3]

	coefficient16 := m[2][0]*m[3][2] - m[3][0]*m[2][2]
	coefficient18 := m[1][0]*m[3][2] - m[3][0]*m[1][2]
	coefficient19 := m[1][0]*m[2][2] - m[2][0]*m[1][2]

	coefficient20 := m[2][0]*m[3][1] - m[3][0]*m[2][1]
	coefficient22 := m[1][0]*m[3][1] - m[3][0]*m[1][1]
	coefficient23 := m[1][0]*m[2][1] - m[2][0]*m[1][1]

	m00 := m[1][0]*coefficient0 - m[1][2]*coefficient4 + m[1][3]*coefficient8
	m01 := -m[0][1]*coefficient0 + m[0][2]*coefficient4 - m[0][3]*coefficient8
	m02 := m[0][1]*coefficient2 - m[0][2]*coefficient6 + m[0][3]*coefficient10
	m03 := -m[0][1]*coefficient3 + m[0][2]*coefficient7 - m[0][3]*coefficient11
	m10 := -m[1][0]*coefficient0 + m[1][2]*coefficient12 - m[1][3]*coefficient16
	m11 := m[0][0]*coefficient0 - m[0][2]*coefficient12 + m[0][3]*coefficient16
	m12 := -m[0][0]*coefficient2 + m[0][2]*coefficient14 - m[0][3]*coefficient18
	m13 := m[0][0]*coefficient3 - m[0][2]*coefficient16 + m[0][3]*coefficient19
	m20 := m[1][0]*coefficient16 - m[1][1]*coefficient12 + m[1][3]*coefficient20
	m21 := -m[0][0]*coefficient16 + m[0][1]*coefficient12 - m[0][3]*coefficient20
	m22 := m[0][0]*coefficient18 - m[0][1]*coefficient14 + m[0][3]*coefficient22
	m23 := -m[0][0]*coefficient19 + m[0][1]*coefficient15 - m[0][3]*coefficient23
	m30 := -m[1][0]*coefficient8 + m[1][1]*coefficient16 - m[1][2]*coefficient20
	m31 := m[0][0]*coefficient8 - m[0][1]*coefficient16 + m[0][2]*coefficient20
	m32 := -m[0][0]*coefficient10 + m[0][1]*coefficient18 - m[0][2]*coefficient22
	m33 := m[0][0]*coefficient11 - m[0][1]*coefficient19 + m[0][2]*coefficient23

	determinant := m[0][0]*m00 + m[0][1]*m10 + m[0][2]*m20 + m[0][3]*m30
	oneOverDeterminant := 1.0 / determinant

	m[0][0] = m00 * oneOverDeterminant
	m[0][1] = m01 * oneOverDeterminant
	m[0][2] = m02 * oneOverDeterminant
	m[0][3] = m03 * oneOverDeterminant
	m[1][0] = m10 * oneOverDeterminant
	m[1][1] = m11 * oneOverDeterminant
	m[1][2] = m12 * oneOverDeterminant
	m[1][3] = m13 * oneOverDeterminant
	m[2][0] = m20 * oneOverDeterminant
	m[2][1] = m21 * oneOverDeterminant
	m[2][2] = m22 * oneOverDeterminant
	m[2][3] = m23 * oneOverDeterminant
	m[3][0] = m30 * oneOverDeterminant
	m[3][1] = m31 * oneOverDeterminant
	m[3][2] = m32 * oneOverDeterminant
	m[3][3] = m33 * oneOverDeterminant

	return m
}

func (m *Mat4x4[T]) MultMatrix4x4(other *Mat4x4[T]) *Mat4x4[T] {
	m00 := m[0][0]*other[0][0] + m[1][0]*other[0][1] + m[2][0]*other[0][2] + m[3][0]*other[0][3]
	m10 := m[0][0]*other[1][0] + m[1][0]*other[1][1] + m[2][0]*other[1][2] + m[3][0]*other[1][3]
	m20 := m[0][0]*other[2][0] + m[1][0]*other[2][1] + m[2][0]*other[2][2] + m[3][0]*other[2][3]
	m30 := m[0][0]*other[3][0] + m[1][0]*other[3][1] + m[2][0]*other[3][2] + m[3][0]*other[3][3]

	m01 := m[0][1]*other[0][0] + m[1][1]*other[0][1] + m[2][1]*other[0][2] + m[3][1]*other[0][3]
	m11 := m[0][1]*other[1][0] + m[1][1]*other[1][1] + m[2][1]*other[1][2] + m[3][1]*other[1][3]
	m21 := m[0][1]*other[2][0] + m[1][1]*other[2][1] + m[2][1]*other[2][2] + m[3][1]*other[2][3]
	m31 := m[0][1]*other[3][0] + m[1][1]*other[3][1] + m[2][1]*other[3][2] + m[3][1]*other[3][3]

	m02 := m[0][2]*other[0][0] + m[1][2]*other[0][1] + m[2][2]*other[0][2] + m[3][2]*other[0][3]
	m12 := m[0][2]*other[1][0] + m[1][2]*other[1][1] + m[2][2]*other[1][2] + m[3][2]*other[1][3]
	m22 := m[0][2]*other[2][0] + m[1][2]*other[2][1] + m[2][2]*other[2][2] + m[3][2]*other[2][3]
	m32 := m[0][2]*other[3][0] + m[1][2]*other[3][1] + m[2][2]*other[3][2] + m[3][2]*other[3][3]

	m03 := m[0][3]*other[0][0] + m[1][3]*other[0][1] + m[2][3]*other[0][2] + m[3][3]*other[0][3]
	m13 := m[0][3]*other[1][0] + m[1][3]*other[1][1] + m[2][3]*other[1][2] + m[3][3]*other[1][3]
	m23 := m[0][3]*other[2][0] + m[1][3]*other[2][1] + m[2][3]*other[2][2] + m[3][3]*other[2][3]
	m33 := m[0][3]*other[3][0] + m[1][3]*other[3][1] + m[2][3]*other[3][2] + m[3][3]*other[3][3]

	m[0][0] = m00
	m[0][1] = m01
	m[0][2] = m02
	m[0][3] = m03
	m[1][0] = m10
	m[1][1] = m11
	m[1][2] = m12
	m[1][3] = m13
	m[2][0] = m20
	m[2][1] = m21
	m[2][2] = m22
	m[2][3] = m23
	m[3][0] = m30
	m[3][1] = m31
	m[3][2] = m32
	m[3][3] = m33

	return m
}

func (m *Mat4x4[T]) Proj3D(normal *Vec3[T]) *Mat4x4[T] {
	var proj Mat4x4[T]
	proj[0][0] = 1 - normal.X*normal.X
	proj[0][1] = -normal.X * normal.Y
	proj[0][2] = -normal.X * normal.Z
	proj[0][3] = 0
	proj[1][0] = -normal.X * normal.Y
	proj[1][1] = 1 - normal.Y*normal.Y
	proj[1][2] = -normal.Y * normal.Z
	proj[1][3] = 0
	proj[2][0] = -normal.X * normal.Z
	proj[2][1] = -normal.Y * normal.Z
	proj[2][2] = 1 - normal.Z*normal.Z
	proj[2][3] = 0
	proj[3][0] = 0
	proj[3][1] = 0
	proj[3][2] = 0
	proj[3][3] = 1

	return m.MultMatrix4x4(&proj)
}

func (m *Mat4x4[T]) Translate(v *Vec3[T]) *Mat4x4[T] {
	m[3][0] = m[0][0]*v.X + m[1][0]*v.Y + m[2][0]*v.Z + m[3][0]
	m[3][1] = m[0][1]*v.X + m[1][1]*v.Y + m[2][1]*v.Z + m[3][1]
	m[3][2] = m[0][2]*v.X + m[1][2]*v.Y + m[2][2]*v.Z + m[3][2]
	m[3][3] = m[0][3]*v.X + m[1][3]*v.Y + m[2][3]*v.Z + m[3][3]

	return m
}

func (m *Mat4x4[T]) Scale(v *Vec3[T]) *Mat4x4[T] {
	m[0][0] *= v.X
	m[0][1] *= v.X
	m[0][2] *= v.X
	m[0][3] *= v.X
	m[1][0] *= v.Y
	m[1][1] *= v.Y
	m[1][2] *= v.Y
	m[1][3] *= v.Y
	m[2][0] *= v.Z
	m[2][1] *= v.Z
	m[2][2] *= v.Z
	m[2][3] *= v.Z

	return m
}

func (m *Mat4x4[T]) RotateAroundAxis(axis *Vec3[T], angleRad T) *Mat4x4[T] {
	var unitAxis Vec3[T]
	unitAxis.SetVec3(axis).Normalize()

	cos := T(math.Cos(float64(angleRad)))
	sin := T(math.Sin(float64(angleRad)))

	inverseCos := 1 - cos

	rotate00 := cos + unitAxis.X*unitAxis.X*inverseCos
	rotate01 := unitAxis.X*unitAxis.Y*inverseCos + sin*unitAxis.Z
	rotate02 := unitAxis.X*unitAxis.Z*inverseCos - sin*unitAxis.Y

	rotate10 := unitAxis.Y*unitAxis.X*inverseCos - sin*unitAxis.Z
	rotate11 := cos + unitAxis.Y*unitAxis.Y*inverseCos
	rotate12 := unitAxis.Y*unitAxis.Z*inverseCos + sin*unitAxis.X

	rotate20 := unitAxis.Z*unitAxis.X*inverseCos + sin*unitAxis.Y
	rotate21 := unitAxis.Z*unitAxis.Y*inverseCos - sin*unitAxis.X
	rotate22 := cos + unitAxis.Z*unitAxis.Z*inverseCos

	m00 := m[0][0]*rotate00 + m[1][0]*rotate01 + m[2][0]*rotate02
	m01 := m[0][1]*rotate00 + m[1][1]*rotate01 + m[2][1]*rotate02
	m02 := m[0][2]*rotate00 + m[1][2]*rotate01 + m[2][2]*rotate02
	m03 := m[0][3]*rotate00 + m[1][3]*rotate01 + m[2][3]*rotate02

	m10 := m[0][0]*rotate10 + m[1][0]*rotate11 + m[2][0]*rotate12
	m11 := m[0][1]*rotate10 + m[1][1]*rotate11 + m[2][1]*rotate12
	m12 := m[0][2]*rotate10 + m[1][2]*rotate11 + m[2][2]*rotate12
	m13 := m[0][3]*rotate10 + m[1][3]*rotate11 + m[2][3]*rotate12

	m20 := m[0][0]*rotate20 + m[1][0]*rotate21 + m[2][0]*rotate22
	m21 := m[0][1]*rotate20 + m[1][1]*rotate21 + m[2][1]*rotate22
	m22 := m[0][2]*rotate20 + m[1][2]*rotate21 + m[2][2]*rotate22
	m23 := m[0][3]*rotate20 + m[1][3]*rotate21 + m[2][3]*rotate22

	m[0][0] = m00
	m[0][1] = m01
	m[0][2] = m02
	m[0][3] = m03
	m[1][0] = m10
	m[1][1] = m11
	m[1][2] = m12
	m[1][3] = m13
	m[2][0] = m20
	m[2][1] = m21
	m[2][2] = m22
	m[2][3] = m23

	return m
}

func (m *Mat4x4[T]) RotateX(angleRad T) *Mat4x4[T] {
	cos := T(math.Cos(float64(angleRad)))
	sin := T(math.Sin(float64(angleRad)))

	m10 := m[1][0]*cos + m[2][0]*sin
	m11 := m[1][1]*cos + m[2][1]*sin
	m12 := m[1][2]*cos + m[2][2]*sin
	m13 := m[1][3]*cos + m[2][3]*sin

	m20 := -m[1][0]*sin + m[2][0]*cos
	m21 := -m[1][1]*sin + m[2][1]*cos
	m22 := -m[1][2]*sin + m[2][2]*cos
	m23 := -m[1][3]*sin + m[2][3]*cos

	m[1][0] = m10
	m[1][1] = m11
	m[1][2] = m12
	m[1][3] = m13
	m[2][0] = m20
	m[2][1] = m21
	m[2][2] = m22
	m[2][3] = m23

	return m
}

func (m *Mat4x4[T]) RotateY(angleRad T) *Mat4x4[T] {
	cos := T(math.Cos(float64(angleRad)))
	sin := T(math.Sin(float64(angleRad)))

	m00 := m[0][0]*cos - m[2][0]*sin
	m01 := m[0][1]*cos - m[2][1]*sin
	m02 := m[0][2]*cos - m[2][2]*sin
	m03 := m[0][3]*cos - m[2][3]*sin

	m20 := m[0][0]*sin + m[2][0]*cos
	m21 := m[0][1]*sin + m[2][1]*cos
	m22 := m[0][2]*sin + m[2][2]*cos
	m23 := m[0][3]*sin + m[2][3]*cos

	m[0][0] = m00
	m[0][1] = m01
	m[0][2] = m02
	m[0][3] = m03
	m[2][0] = m20
	m[2][1] = m21
	m[2][2] = m22
	m[2][3] = m23

	return m
}

func (m *Mat4x4[T]) RotateZ(angleRad T) *Mat4x4[T] {
	cos := T(math.Cos(float64(angleRad)))
	sin := T(math.Sin(float64(angleRad)))

	m00 := m[0][0]*cos + m[1][0]*sin
	m01 := m[0][1]*cos + m[1][1]*sin
	m02 := m[0][2]*cos + m[1][2]*sin
	m03 := m[0][3]*cos + m[1][3]*sin

	m10 := -m[0][0]*sin + m[1][0]*cos
	m11 := -m[0][1]*sin + m[1][1]*cos
	m12 := -m[0][2]*sin + m[1][2]*cos
	m13 := -m[0][3]*sin + m[1][3]*cos

	m[0][0] = m00
	m[0][1] = m01
	m[0][2] = m02
	m[0][3] = m03
	m[1][0] = m10
	m[1][1] = m11
	m[1][2] = m12
	m[1][3] = m13

	return m
}

func (m *Mat4x4[T]) ShearX(y, z T) *Mat4x4[T] {

	m01 := m[0][0]*y + m[0][1]
	m02 := m[0][0]*z + m[0][2]

	m11 := m[1][0]*y + m[1][1]
	m12 := m[1][0]*z + m[1][2]

	m21 := m[2][0]*y + m[2][1]
	m22 := m[2][0]*z + m[2][2]

	m31 := m[3][0]*y + m[3][1]
	m32 := m[3][0]*z + m[3][2]

	m[0][1] = m01
	m[0][2] = m02
	m[1][1] = m11
	m[1][2] = m12
	m[2][1] = m21
	m[2][2] = m22
	m[3][1] = m31
	m[3][2] = m32

	return m
}

func (m *Mat4x4[T]) ShearY(x, z T) *Mat4x4[T] {
	m00 := m[0][1]*x + m[0][0]
	m02 := m[0][1]*z + m[0][2]

	m10 := m[1][1]*x + m[1][0]
	m12 := m[1][1]*z + m[1][2]

	m20 := m[2][1]*x + m[2][0]
	m22 := m[2][1]*z + m[2][2]

	m[0][0] = m00
	m[0][2] = m02
	m[1][0] = m10
	m[1][2] = m12
	m[2][0] = m20
	m[2][2] = m22

	return m
}

func (m *Mat4x4[T]) ShearZ(x, y T) *Mat4x4[T] {
	m00 := m[0][2]*x + m[0][0]
	m01 := m[0][2]*y + m[0][1]

	m10 := m[1][2]*x + m[1][0]
	m11 := m[1][2]*y + m[1][1]

	m20 := m[2][2]*x + m[2][0]
	m21 := m[2][2]*y + m[2][1]

	m[0][0] = m00
	m[0][1] = m01
	m[1][0] = m10
	m[1][1] = m11
	m[2][0] = m20
	m[2][1] = m21

	return m
}

func (m *Mat4x4[T]) InterpolateMatrix(otherMatrix *Mat4x4[T], delta T) *Mat4x4[T] {
	var oldMatrix Mat4x4[T]
	oldMatrix.SetMat4x4(m)

	var thisRotation Mat4x4[T]
	thisRotation.SetMatrixRotationFrom(m)

	var transposedRotation Mat4x4[T]
	transposedRotation.SetMat4x4(&thisRotation).Transpose()

	var deltaRotation Mat4x4[T]
	deltaRotation.SetMultMatrix4x4(otherMatrix, &transposedRotation)

	var deltaAxis Vec3[T]
	var deltaAngle T
	deltaRotation.GetAxisAngle(&deltaAxis, &deltaAngle)

	m.SetAxisAngle(&deltaAxis, deltaAngle).MultMatrix4x4(&thisRotation)
	m[3][0] = oldMatrix[3][0] + delta*(otherMatrix[3][0]-oldMatrix[3][0])
	m[3][1] = oldMatrix[3][1] + delta*(otherMatrix[3][1]-oldMatrix[3][1])
	m[3][2] = oldMatrix[3][2] + delta*(otherMatrix[3][2]-oldMatrix[3][2])

	return m
}

func (m *Mat4x4[T]) IsNormalized(epsilon T) bool {
	for i := 0; i < 4; i++ {
		column := abs[T](m[i][0]*m[i][0] + m[i][1]*m[i][1] + m[i][2]*m[i][2] + m[i][3]*m[i][3])
		if column-T(1) > T(2)*epsilon {
			return false
		}
	}

	for i := 0; i < 4; i++ {
		row := abs[T](m[0][i]*m[0][i] + m[1][i]*m[1][i] + m[2][i]*m[2][i] + m[3][i]*m[3][i])
		for row-T(1) > T(2)*epsilon {
			return false
		}
	}

	return true
}

func (m *Mat4x4[T]) IsNull(epsilon T) bool {
	for i := 0; i < 4; i++ {
		column := abs[T](m[i][0]*m[i][0] + m[i][1]*m[i][1] + m[i][2]*m[i][2] + m[i][3]*m[i][3])
		if column > epsilon {
			return false
		}
	}

	return true
}

func (m *Mat4x4[T]) Equal(other *Mat4x4[T], epsilon T) bool {
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			if abs[T](m[i][j]-other[i][j]) > epsilon {
				return false
			}
		}
	}

	return true
}
