package math

import "math"

type Mat4Row [4]float32
type Mat4x4 [4]Mat4Row

func (m *Mat4x4) SetIdentity() *Mat4x4 {
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

func (m *Mat4x4) SetDiagonalScalar(s float32) *Mat4x4 {
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

func (m *Mat4x4) SetDiagonalVector(v *Vec4) *Mat4x4 {
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

func (m *Mat4x4) SetFrustum(left, right, bottom, top, nearVal, farVal float32) *Mat4x4 {
	m[0][0] = (2.0 * nearVal) / (right - left)
	m[0][1] = 0
	m[0][2] = 0
	m[0][3] = 0
	m[1][0] = 0
	m[1][1] = (2 * nearVal) / (top - bottom)
	m[1][2] = 0
	m[1][3] = 0
	m[2][0] = (right + left) / (right - left)
	m[2][1] = (top + bottom) / (top - bottom)
	m[2][2] = farVal / (nearVal - farVal)
	m[2][3] = -1
	m[3][0] = 0
	m[3][1] = 0
	m[3][2] = -(farVal * nearVal) / (farVal - nearVal)
	m[3][3] = 0

	return m
}

func (m *Mat4x4) SetAxisAngle(axis *Vec3, angleRad float32) *Mat4x4 {
	cos := float32(math.Cos(float64(angleRad)))
	sin := float32(math.Sin(float64(angleRad)))
	inverseCos := 1 - cos

	var unitAxis Vec3
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

func (m *Mat4x4) SetOrientation(normal *Vec3, up *Vec3) *Mat4x4 {
	if normal.Equal(up, 0.0001) {
		m.SetIdentity()
		return m
	}

	var rotationAxis Vec3
	rotationAxis.SetVec3(up).CrossProduct(normal)

	angle := float32(math.Acos(float64(normal.DotProduct(up))))

	return m.SetAxisAngle(&rotationAxis, angle)
}

func (m *Mat4x4) SetTranslation(x, y, z float32) *Mat4x4 {
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

func (m *Mat4x4) SetScale(x, y, z float32) *Mat4x4 {
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

func (m *Mat4x4) SetMat4x4(other *Mat4x4) *Mat4x4 {
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

func (m *Mat4x4) SetMat3x3(other *Mat3x3) *Mat4x4 {
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

func (m *Mat4x4) SetColMajor(c1, c2, c3, c4 *Vec4) *Mat4x4 {
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

func (m *Mat4x4) SetRowMajor(r1, r2, r3, r4 *Vec4) *Mat4x4 {
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

func (m *Mat4x4) SetEulerAngles(yawRad, pitchRad, rollRad float32) *Mat4x4 {
	yawCos := float32(math.Cos(float64(yawRad)))
	pitchCos := float32(math.Cos(float64(pitchRad)))
	rollCos := float32(math.Cos(float64(rollRad)))
	yawSin := float32(math.Sin(float64(yawRad)))
	pitchSin := float32(math.Sin(float64(pitchRad)))
	rollSin := float32(math.Sin(float64(rollRad)))

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

func (m *Mat4x4) SetMatrixRotationFrom(other *Mat4x4) *Mat4x4 {
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

func (m *Mat4x4) SetInfinitePerspective(fovY, aspectRatio, zNear float32) *Mat4x4 {
	r := float32(math.Tan(float64(fovY/2))) * zNear
	left := -r * aspectRatio
	right := r * aspectRatio
	bottom := -r
	top := r

	m[0][0] = (2 * zNear) / (right - left)
	m[0][1] = 0
	m[0][2] = 0
	m[0][3] = 0
	m[1][0] = 0
	m[1][1] = (2 * zNear) / (top - bottom)
	m[1][2] = 0
	m[1][3] = 0
	m[2][0] = 0
	m[2][1] = 0
	m[2][2] = -1
	m[2][3] = -1
	m[3][0] = 0
	m[3][1] = 0
	m[3][2] = -2 * zNear
	m[3][3] = 0

	return m
}

func (m *Mat4x4) SetOrthographic2D(left, right, bottom, top float32) *Mat4x4 {
	m[0][0] = 2 / (right - left)
	m[0][1] = 0
	m[0][2] = 0
	m[0][3] = 0
	m[1][0] = 0
	m[1][1] = 2 / (top - bottom)
	m[1][2] = 0
	m[1][3] = 0
	m[2][0] = 0
	m[2][1] = 0
	m[2][2] = 1
	m[2][3] = 0
	m[3][0] = -(right + left) / (right - left)
	m[3][1] = -(top + bottom) / (top - bottom)
	m[3][2] = 0
	m[3][3] = 1

	return m
}

func (m *Mat4x4) SetOrthographic(left, right, bottom, top, zNear, zFar float32) *Mat4x4 {
	m[0][0] = 2 / (right - left)
	m[0][1] = 0
	m[0][2] = 0
	m[0][3] = 0
	m[1][0] = 0
	m[1][1] = 2 / (top - bottom)
	m[1][2] = 0
	m[1][3] = 0
	m[2][0] = 0
	m[2][1] = 0
	m[2][2] = -1 / (zFar - zNear)
	m[2][3] = 0
	m[3][0] = -(right + left) / (right - left)
	m[3][1] = -(top + bottom) / (top - bottom)
	m[3][2] = -zNear / (zFar - zNear)
	m[3][3] = 1

	return m
}

func (m *Mat4x4) SetPerspective(fovY, aspectRatio, zNear, zFar float32) *Mat4x4 {
	tanHalfFovY := float32(math.Tan(float64(fovY / 2)))

	m[0][0] = 1 / (aspectRatio * tanHalfFovY)
	m[0][1] = 0
	m[0][2] = 0
	m[0][3] = 0
	m[1][0] = 0
	m[1][1] = 1 / tanHalfFovY
	m[1][2] = 0
	m[1][3] = 0
	m[2][0] = 0
	m[2][1] = 0
	m[2][2] = zFar / (zNear - zFar)
	m[2][3] = -1
	m[3][0] = 0
	m[3][1] = 0
	m[3][2] = -(zFar * zNear) / (zFar - zNear)
	m[3][3] = 0

	return m
}

func (m *Mat4x4) SetPerspectiveFOV(fov, width, height, zNear, zFar float32) *Mat4x4 {
	h := float32(math.Cos(0.5*float64(fov))) / float32(math.Sin(0.5*float64(fov)))
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
	m[2][3] = -1
	m[3][0] = 0
	m[3][1] = 0
	m[3][2] = -(zFar * zNear) / (zFar - zNear)
	m[3][3] = 0

	return m
}

func (m *Mat4x4) SetPickMatrix(center *Vec2, delta *Vec2, viewport Vec4) *Mat4x4 {
	m.SetIdentity()

	if !(delta.X > 0 && delta.Y > 0) {
		return m
	}

	var temp Vec3
	temp.X = (viewport.Z - 2*(center.X-viewport.X)) / delta.X
	temp.Y = (viewport.W - 2*(center.Y-viewport.Y)) / delta.Y
	temp.Z = 0

	var scale Vec3
	scale.X = viewport.Z / delta.X
	scale.Y = viewport.W / delta.Y
	scale.Z = 1

	return m.Translate(&temp).Scale(&scale)
}

func (m *Mat4x4) SetLookAt(eyePosition *Vec3, target *Vec3, up *Vec3) *Mat4x4 {
	var f Vec3
	f.SetVec3(target).SubtractVec3(eyePosition).Normalize()

	var s Vec3
	s.SetVec3(&f).CrossProduct(up).Normalize()

	var u Vec3
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

func (m *Mat4x4) SetMatrixCrossProduct(vec *Vec3) *Mat4x4 {
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

func (m *Mat4x4) SetMultMatrix4x4(lhs, rhs *Mat4x4) *Mat4x4 {
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

func (m *Mat4x4) SetInterpolate(lhs, rhs *Mat4x4, delta float32) *Mat4x4 {
	return m.SetMat4x4(lhs).InterpolateMatrix(rhs, delta)
}

func (m *Mat4x4) GetAxisAngle(outAxis *Vec3, outAngleRad *float32) {
	epsilon := float32(0.01)
	epsilon2 := float32(0.1)

	if (abs(m[1][0]-m[0][1]) < epsilon) &&
		(abs(m[2][0]-m[0][2]) < epsilon) &&
		(abs(m[2][1]-m[1][2]) < epsilon) {

		if (abs(m[1][0]+m[0][1]) < epsilon2) &&
			(abs(m[2][0]+m[0][2]) < epsilon2) &&
			(abs(m[2][1]+m[1][2]) < epsilon2) &&
			(abs(m[0][0]+m[1][1]+m[2][2]-3) < epsilon2) {

			*outAngleRad = 0
			outAxis.X = 1
			outAxis.Y = 0
			outAxis.Z = 0
			return
		}

		*outAngleRad = float32(math.Pi)
		xx := (m[0][0] + 1) * 0.5
		yy := (m[1][1] + 1) * 0.5
		zz := (m[2][2] + 1) * 0.5
		xy := (m[1][0] + m[0][1]) * 0.25
		xz := (m[2][0] + m[0][2]) * 0.25
		yz := (m[2][1] + m[1][2]) * 0.25

		if xx > yy && xx > zz {
			if xx < epsilon {
				outAxis.X = 0
				outAxis.Y = 0.7071
				outAxis.Z = 0.7071
			} else {
				outAxis.X = float32(math.Sqrt(float64(xx)))
				outAxis.Y = xy / outAxis.X
				outAxis.Z = xz / outAxis.X
			}
		} else if yy > zz {
			if yy < epsilon {
				outAxis.X = 0.7071
				outAxis.Y = 0
				outAxis.Z = 0.7071
			} else {
				outAxis.Y = float32(math.Sqrt(float64(yy)))
				outAxis.X = xy / outAxis.Y
				outAxis.Z = yz / outAxis.Z
			}
		} else {
			if zz < epsilon {
				outAxis.X = 0.7071
				outAxis.Y = 0.7071
				outAxis.Z = 0
			} else {
				outAxis.Z = float32(math.Sqrt(float64(zz)))
				outAxis.X = xz / outAxis.Z
				outAxis.Y = yz / outAxis.Z
			}
		}
		return
	}

	sSquared := (m[2][1]-m[1][2])*(m[2][1]-m[1][2]) +
		(m[2][0]-m[0][2])*(m[2][0]-m[0][2]) +
		(m[1][0]-m[0][1])*(m[1][0]-m[0][1])
	s := float32(math.Sqrt(float64(sSquared)))

	if s < 0.001 {
		s = 1
	}

	angleCos := (m[0][0] + m[1][1] + m[2][2] - 1) * 0.5
	if angleCos-1 < epsilon {
		*outAngleRad = float32(math.Pi) * 0.25
	} else {
		*outAngleRad = float32(math.Acos(float64(angleCos)))
	}

	outAxis.X = (m[1][2] - m[2][1]) / s
	outAxis.Y = (m[2][0] - m[0][2]) / s
	outAxis.Z = (m[0][1] - m[1][0]) / s
}

func (m *Mat4x4) Transpose() *Mat4x4 {
	m[1][0], m[0][1] = m[0][1], m[1][0]
	m[2][0], m[0][2] = m[0][2], m[2][0]
	m[2][1], m[1][2] = m[1][2], m[2][1]

	m[3][0], m[0][3] = m[0][3], m[3][0]
	m[3][1], m[1][3] = m[1][3], m[3][1]
	m[3][2], m[2][3] = m[2][3], m[3][2]

	return m
}

func (m *Mat4x4) Inverse() *Mat4x4 {
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

func (m *Mat4x4) MultMatrix4x4(other *Mat4x4) *Mat4x4 {
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

func (m *Mat4x4) Proj3D(normal *Vec3) *Mat4x4 {
	var proj Mat4x4
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

func (m *Mat4x4) Translate(v *Vec3) *Mat4x4 {
	m[3][0] = m[0][0]*v.X + m[1][0]*v.Y + m[2][0]*v.Z + m[3][0]
	m[3][1] = m[0][1]*v.X + m[1][1]*v.Y + m[2][1]*v.Z + m[3][1]
	m[3][2] = m[0][2]*v.X + m[1][2]*v.Y + m[2][2]*v.Z + m[3][2]
	m[3][3] = m[0][3]*v.X + m[1][3]*v.Y + m[2][3]*v.Z + m[3][3]

	return m
}

func (m *Mat4x4) Scale(v *Vec3) *Mat4x4 {
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

func (m *Mat4x4) RotateAroundAxis(axis *Vec3, angleRad float32) *Mat4x4 {
	var unitAxis Vec3
	unitAxis.SetVec3(axis).Normalize()

	cos := float32(math.Cos(float64(angleRad)))
	sin := float32(math.Sin(float64(angleRad)))

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

func (m *Mat4x4) RotateX(angleRad float32) *Mat4x4 {
	cos := float32(math.Cos(float64(angleRad)))
	sin := float32(math.Sin(float64(angleRad)))

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

func (m *Mat4x4) RotateZ(angleRad float32) *Mat4x4 {
	cos := float32(math.Cos(float64(angleRad)))
	sin := float32(math.Sin(float64(angleRad)))

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

func (m *Mat4x4) ShearX(y, z float32) *Mat4x4 {

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

func (m *Mat4x4) ShearY(x, z float32) *Mat4x4 {
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

func (m *Mat4x4) ShearZ(x, y float32) *Mat4x4 {
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

func (m *Mat4x4) InterpolateMatrix(otherMatrix *Mat4x4, delta float32) *Mat4x4 {
	var oldMatrix Mat4x4
	oldMatrix.SetMat4x4(m)

	var thisRotation Mat4x4
	thisRotation.SetMatrixRotationFrom(m)

	var transposedRotation Mat4x4
	transposedRotation.SetMat4x4(&thisRotation).Transpose()

	var deltaRotation Mat4x4
	deltaRotation.SetMat4x4(otherMatrix).MultMatrix4x4(&transposedRotation)

	var deltaAxis Vec3
	var deltaAngle float32
	deltaRotation.GetAxisAngle(&deltaAxis, &deltaAngle)

	m.SetAxisAngle(&deltaAxis, deltaAngle).MultMatrix4x4(&thisRotation)
	m[3][0] = oldMatrix[3][0] + delta*(otherMatrix[3][0]-oldMatrix[3][0])
	m[3][1] = oldMatrix[3][1] + delta*(otherMatrix[3][1]-oldMatrix[3][1])
	m[3][2] = oldMatrix[3][2] + delta*(otherMatrix[3][2]-oldMatrix[3][2])

	return m
}

func (m *Mat4x4) IsNormalized(epsilon float32) bool {
	for i := 0; i < 4; i++ {
		column := abs(m[i][0]*m[i][0] + m[i][1]*m[i][1] + m[i][2]*m[i][2] + m[i][3]*m[i][3])
		if column-1 > 2*epsilon {
			return false
		}
	}

	for i := 0; i < 4; i++ {
		row := abs(m[0][i]*m[0][i] + m[1][i]*m[1][i] + m[2][i]*m[2][i] + m[3][i]*m[3][i])
		for row-1 > 2*epsilon {
			return false
		}
	}

	return true
}

func (m *Mat4x4) IsNull(epsilon float32) bool {
	for i := 0; i < 4; i++ {
		column := abs(m[i][0]*m[i][0] + m[i][1]*m[i][1] + m[i][2]*m[i][2] + m[i][3]*m[i][3])
		if column > epsilon {
			return false
		}
	}

	return true
}

func (m *Mat4x4) Equal(other *Mat4x4, epsilon float32) bool {
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			if abs(m[i][j]-other[i][j]) > epsilon {
				return false
			}
		}
	}

	return true
}
