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

func (m *Mat4x4[T]) SetDiagonal(v *Vec4[T]) *Mat4x4[T] {
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
	t := T(1) - cos

	var unitAxis Vec3[T]
	unitAxis.CopyFrom(axis).Normalize()

	m[0][0] = t*unitAxis.X*unitAxis.X + cos
	m[0][1] = t*unitAxis.X*unitAxis.Y + unitAxis.Z*sin
	m[0][2] = t*unitAxis.X*unitAxis.Z - unitAxis.Y*sin
	m[0][3] = 0
	m[1][0] = t*unitAxis.X*unitAxis.Y - unitAxis.Z*sin
	m[1][1] = t*unitAxis.Y*unitAxis.Y + cos
	m[1][2] = t*unitAxis.Y*unitAxis.Z + unitAxis.X*sin
	m[1][3] = 0
	m[2][0] = t*unitAxis.X*unitAxis.Z + unitAxis.Y*sin
	m[2][1] = t*unitAxis.Y*unitAxis.Z - unitAxis.X*sin
	m[2][2] = t*unitAxis.Z*unitAxis.Z + cos
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
	rotationAxis.CopyFrom(up).CrossProduct(normal)

	angle := T(math.Acos(float64(normal.DotProduct(up))))

	return m.SetAxisAngle(&rotationAxis, angle)
}

func (m *Mat4x4[T]) CopyFrom(other *Mat4x4[T]) *Mat4x4[T] {
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

func (m *Mat4x4[T]) SetEulerAngleYXZ(yawRad, pitchRad, rollRad T) *Mat4x4[T] {
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
	f.CopyFrom(target).SubtractVec3(eyePosition).Normalize()

	var s Vec3[T]
	s.CopyFrom(&f).CrossProduct(up).Normalize()

	var u Vec3[T]
	u.CopyFrom(&s).CrossProduct(&f)

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

func abs[T FloatingPoint](in T) T {
	return T(math.Abs(float64(in)))
}

func (m *Mat4x4[T]) AxisAngle(outAxis *Vec3[T], outAngleRad *T) {
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

func (m *Mat4x4[T]) MultMatrix(other *Mat4x4[T]) *Mat4x4[T] {
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

	return m.MultMatrix(&proj)
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

func (m *Mat4x4[T]) RotateAroundPrenormalizedAxis(unitAxis *Vec3[T], angleRad T) *Mat4x4[T] {
	cos := T(math.Cos(float64(angleRad)))
	sin := T(math.Sin(float64(angleRad)))

	var temp Vec3[T]
	temp.CopyFrom(unitAxis).Scale(1 - cos)

	rotate00 := cos + temp.X*unitAxis.X
	rotate01 := temp.X*unitAxis.Y + sin*unitAxis.Z
	rotate02 := temp.X*unitAxis.Z - sin*unitAxis.Y

	rotate10 := temp.Y*unitAxis.X - sin*unitAxis.Z
	rotate11 := cos + temp.Y*unitAxis.Y
	rotate12 := temp.Y*unitAxis.Z + sin*unitAxis.X

	rotate20 := temp.Z*unitAxis.X + sin*unitAxis.Y
	rotate21 := temp.Z*unitAxis.Y - sin*unitAxis.X
	rotate22 := cos + temp.Z*unitAxis.Z

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

func (m *Mat4x4[T]) RotateAroundAxis(axis *Vec3[T], angleRad T) *Mat4x4[T] {
	var unitAxis Vec3[T]
	unitAxis.CopyFrom(axis).Normalize()

	return m.RotateAroundPrenormalizedAxis(&unitAxis, angleRad)
}

func (m *Mat4x4[T]) RotateX(angleRad T) *Mat4x4[T] {
	axis := Vec3[T]{X: 1, Y: 0, Z: 0}
	return m.RotateAroundPrenormalizedAxis(&axis, angleRad)
}

func (m *Mat4x4[T]) RotateY(angleRad T) *Mat4x4[T] {
	axis := Vec3[T]{X: 0, Y: 1, Z: 0}
	return m.RotateAroundPrenormalizedAxis(&axis, angleRad)
}

func (m *Mat4x4[T]) RotateZ(angleRad T) *Mat4x4[T] {
	axis := Vec3[T]{X: 0, Y: 0, Z: 1}
	return m.RotateAroundPrenormalizedAxis(&axis, angleRad)
}

func (m *Mat4x4[T]) InterpolateMatrix(otherMatrix *Mat4x4[T], delta T) *Mat4x4[T] {
	var oldMatrix Mat4x4[T]
	oldMatrix.CopyFrom(m)

	var thisRotation Mat4x4[T]
	thisRotation.SetMatrixRotationFrom(m)

	var transposedRotation Mat4x4[T]
	transposedRotation.CopyFrom(&thisRotation).Transpose()

	var deltaRotation Mat4x4[T]
	deltaRotation.CopyFrom(otherMatrix).MultMatrix(&transposedRotation)

	var deltaAxis Vec3[T]
	var deltaAngle T
	deltaRotation.AxisAngle(&deltaAxis, &deltaAngle)

	m.SetAxisAngle(&deltaAxis, deltaAngle).MultMatrix(&thisRotation)
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
