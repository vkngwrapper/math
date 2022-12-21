package math

import "math"

// Quaternion is a 4-element vector that represents a 3d rotation in a way that is immune
// to gimbal lock
type Quaternion[T FloatingPoint] struct {
	X, Y, Z, W T
}

// SetIdentity overwrites the contents of this quaternion with the identity quaternion
func (q *Quaternion[T]) SetIdentity() *Quaternion[T] {
	q.X = 0
	q.Y = 0
	q.Z = 0
	q.W = 1

	return q
}

// SetRotationAroundAxis overwrites the current contents of this quaternion with a quaternion that rotates
// around the provided axis by the provided angle in radians.
//
// axis - A 3-element vector that is normal to the angle of rotation. It does not need to be normalized.
//
// angleRad - The amount to rotate in radians
func (q *Quaternion[T]) SetRotationAroundAxis(axis *Vec3[T], angle T) *Quaternion[T] {
	var unitAxis Vec3[T]
	unitAxis.SetVec3(axis).Normalize()

	s := T(math.Sin(float64(angle * 0.5)))
	q.X = unitAxis.X * s
	q.Y = unitAxis.Y * s
	q.Z = unitAxis.Z * s
	q.W = T(math.Cos(float64(angle * 0.5)))

	return q
}

// SetRotationX overwrites the current contents of this quaternion with a quaternion that
// rotates around the x axis by the specified amount
//
// pitchRad - The angle to rotate around the x axis in radians
func (q *Quaternion[T]) SetRotationX(angle T) *Quaternion[T] {
	q.X = T(math.Sin(float64(angle * 0.5)))
	q.Y = 0
	q.Z = 0
	q.W = T(math.Cos(float64(angle * 0.5)))

	return q
}

// SetRotationY overwrites the current contents of this quaternion with a quaternion that
// rotates around the y axis by the specified amount
//
// yawRad - The angle to rotate around the y axis in radians
func (q *Quaternion[T]) SetRotationY(angle T) *Quaternion[T] {
	q.X = 0
	q.Y = T(math.Sin(float64(angle * 0.5)))
	q.Z = 0
	q.W = T(math.Cos(float64(angle * 0.5)))

	return q
}

// SetRotationZ overwrites the current contents of this quaternion with a quaternion that
// rotates around the z axis by the specified amount
//
// rollRad - The angle to rotate around the z axis in radians
func (q *Quaternion[T]) SetRotationZ(angle T) *Quaternion[T] {
	q.X = 0
	q.Y = 0
	q.Z = T(math.Sin(float64(angle * 0.5)))
	q.W = T(math.Cos(float64(angle * 0.5)))

	return q
}

// SetOrientation overwrites the current quaternion with a quaternion which rotates
// from the provided source vector to the provided target vector
//
// origin - A unit vector indicating the starting direction of the rotation
//
// target - A unit vector indicating the target direction of the rotation
func (q *Quaternion[T]) SetOrientation(origin *Vec3[T], target *Vec3[T]) *Quaternion[T] {
	cosTheta := origin.DotProduct(target)

	if cosTheta >= 0.9999 {
		q.X = 0
		q.Y = 0
		q.Z = 0
		q.W = 1
		return q
	}

	if cosTheta < -0.9999 {
		rotationAxis := Vec3[T]{0, 0, 1}
		rotationAxis.CrossProduct(origin)
		if rotationAxis.LenSqr() < 0.0001 {
			rotationAxis = Vec3[T]{1, 0, 0}
			rotationAxis.CrossProduct(origin)
		}

		rotationAxis.Normalize()
		q.X = rotationAxis.X
		q.Y = rotationAxis.Y
		q.Z = rotationAxis.Z
		q.W = T(math.Cos(math.Pi * 0.5))
		return q
	}

	var rotationAxis Vec3[T]
	rotationAxis.SetCrossProduct(origin, target)

	lengthSquared := (1.0 + cosTheta) * 2.0
	length := T(math.Sqrt(float64(lengthSquared)))
	inverse := 1.0 / length

	q.X = rotationAxis.X * inverse
	q.Y = rotationAxis.Y * inverse
	q.Z = rotationAxis.Z * inverse
	q.W = length * 0.5

	return q
}

// SetQuaternion overwrites the current quaternion with the contents of the provided quaternion
//
// other - The quaternion to initialize from
func (q *Quaternion[T]) SetQuaternion(other *Quaternion[T]) *Quaternion[T] {
	q.X = other.X
	q.Y = other.Y
	q.Z = other.Z
	q.W = other.W

	return q
}

// SetMat3x3 overwrites the current quaternion with the rotation of a provided transform matrix
//
// m - The transform matrix to initialize from
func (q *Quaternion[T]) SetMat3x3(m *Mat3x3[T]) *Quaternion[T] {
	fourXSquaredMinus1 := m[0][0] - m[1][1] - m[2][2]
	fourYSquaredMinus1 := m[1][1] - m[0][0] - m[2][2]
	fourZSquaredMinus1 := m[2][2] - m[0][0] - m[1][1]
	fourWSquaredMinus1 := m[0][0] + m[1][1] + m[2][2]

	biggestIndex := 0
	fourBiggestSquaredMinus1 := fourWSquaredMinus1
	if fourXSquaredMinus1 > fourBiggestSquaredMinus1 {
		fourBiggestSquaredMinus1 = fourXSquaredMinus1
		biggestIndex = 1
	}
	if fourYSquaredMinus1 > fourBiggestSquaredMinus1 {
		fourBiggestSquaredMinus1 = fourYSquaredMinus1
		biggestIndex = 2
	}
	if fourZSquaredMinus1 > fourBiggestSquaredMinus1 {
		fourBiggestSquaredMinus1 = fourZSquaredMinus1
		biggestIndex = 3
	}

	biggestVal := T(math.Sqrt(float64(fourBiggestSquaredMinus1+1.0))) * 0.5
	mult := 0.25 / biggestVal

	switch biggestIndex {
	case 0:
		q.X = (m[1][2] - m[2][1]) * mult
		q.Y = (m[2][0] - m[0][2]) * mult
		q.Z = (m[0][1] - m[1][0]) * mult
		q.W = biggestVal
		return q
	case 1:
		q.X = biggestVal
		q.Y = (m[0][1] + m[1][0]) * mult
		q.Z = (m[2][0] + m[0][2]) * mult
		q.W = (m[1][2] - m[2][1]) * mult
		return q
	case 2:
		q.X = (m[0][1] + m[1][0]) * mult
		q.Y = biggestVal
		q.Z = (m[1][2] + m[2][1]) * mult
		q.W = (m[2][0] - m[0][2]) * mult
		return q
	case 3:
		q.X = (m[2][0] + m[0][2]) * mult
		q.Y = (m[1][2] + m[2][1]) * mult
		q.Z = biggestVal
		q.W = (m[0][1] - m[1][0]) * mult
		return q
	default:
		q.X = 0
		q.Y = 0
		q.Z = 0
		q.W = 1
		return q
	}
}

// SetMat4x4 overwrites the current quaternion with the rotation of a provided transform matrix
//
// m - The transform matrix to initialize from
func (q *Quaternion[T]) SetMat4x4(m *Mat4x4[T]) *Quaternion[T] {
	fourXSquaredMinus1 := m[0][0] - m[1][1] - m[2][2]
	fourYSquaredMinus1 := m[1][1] - m[0][0] - m[2][2]
	fourZSquaredMinus1 := m[2][2] - m[0][0] - m[1][1]
	fourWSquaredMinus1 := m[0][0] + m[1][1] + m[2][2]

	biggestIndex := 0
	fourBiggestSquaredMinus1 := fourWSquaredMinus1
	if fourXSquaredMinus1 > fourBiggestSquaredMinus1 {
		fourBiggestSquaredMinus1 = fourXSquaredMinus1
		biggestIndex = 1
	}
	if fourYSquaredMinus1 > fourBiggestSquaredMinus1 {
		fourBiggestSquaredMinus1 = fourYSquaredMinus1
		biggestIndex = 2
	}
	if fourZSquaredMinus1 > fourBiggestSquaredMinus1 {
		fourBiggestSquaredMinus1 = fourZSquaredMinus1
		biggestIndex = 3
	}

	biggestVal := T(math.Sqrt(float64(fourBiggestSquaredMinus1+1.0))) * 0.5
	mult := 0.25 / biggestVal

	switch biggestIndex {
	case 0:
		q.X = (m[1][2] - m[2][1]) * mult
		q.Y = (m[2][0] - m[0][2]) * mult
		q.Z = (m[0][1] - m[1][0]) * mult
		q.W = biggestVal
		return q
	case 1:
		q.X = biggestVal
		q.Y = (m[0][1] + m[1][0]) * mult
		q.Z = (m[2][0] + m[0][2]) * mult
		q.W = (m[1][2] - m[2][1]) * mult
		return q
	case 2:
		q.X = (m[0][1] + m[1][0]) * mult
		q.Y = biggestVal
		q.Z = (m[1][2] + m[2][1]) * mult
		q.W = (m[2][0] - m[0][2]) * mult
		return q
	case 3:
		q.X = (m[2][0] + m[0][2]) * mult
		q.Y = (m[1][2] + m[2][1]) * mult
		q.Z = biggestVal
		q.W = (m[0][1] - m[1][0]) * mult
		return q
	default:
		q.X = 0
		q.Y = 0
		q.Z = 0
		q.W = 1
		return q
	}
}

// SetIntermediate overwrites the current quaternion with a SQUAD control point quaternion for
// the control point corresponding to the 'current' rotation keyframe. This control point is
// used as an input for the SetSquad method.
//
// prev - The keyframe before 'current' - this may be the same as 'current' when 'current' is the
// first keyframe in a path
//
// current - The keyframe we are attempting to retrieve the control point for
//
// next - The keyframe after 'current' - this may be the same as 'current' when 'current' is the last
// keyframe in a path
func (q *Quaternion[T]) SetIntermediate(prev *Quaternion[T], current *Quaternion[T], next *Quaternion[T]) *Quaternion[T] {
	var inverseCurrent Quaternion[T]
	inverseCurrent.SetQuaternion(current).Inverse()

	var nextDiff Quaternion[T]
	nextDiff.SetMultQuaternion(next, &inverseCurrent).Log()

	var prevDiff Quaternion[T]
	prevDiff.SetMultQuaternion(prev, &inverseCurrent).Log()

	factor := T(-1.0 / 4.0)
	q.X = (nextDiff.X + prevDiff.X) * factor
	q.Y = (nextDiff.Y + prevDiff.Y) * factor
	q.Z = (nextDiff.Z + prevDiff.Z) * factor
	q.W = (nextDiff.W + prevDiff.W) * factor

	q.Exp().MultQuaternion(current)
	return q
}

// SetMix performs an oriented spherical linear interpolation at constant speed between two quaternions
// and overwrites the current quaternion with the result.
//
// lhs - The origin quaternion in the interpolation
//
// rhs - The target quaternion in the interpolation
//
// delta - The mix factor between the two quaternions- this can be any value: between 0 and 1,
// the interpolation will move linearly from the lhs to the rhs quaternion. Beyond those bounds,
// the interpolation oscillates between the two values
func (q *Quaternion[T]) SetMix(lhs, rhs *Quaternion[T], delta T) *Quaternion[T] {
	cosTheta := lhs.DotProduct(rhs)

	// Perform a linear interpolation when cosTheta is close to 1 to avoid side effect of sin(angle) becoming a zero denominator
	if cosTheta > 0.9999 {
		q.X = mix[T](lhs.X, rhs.X, delta)
		q.Y = mix[T](lhs.Y, rhs.Y, delta)
		q.Z = mix[T](lhs.Z, rhs.Z, delta)
		q.W = mix[T](lhs.W, rhs.W, delta)

		return q
	}

	angle := math.Acos(float64(cosTheta))
	qFactor := T(math.Sin(float64(1.0-delta) * angle))
	otherFactor := T(math.Sin(float64(delta) * angle))
	inverseFactor := 1.0 / T(math.Sin(float64(angle)))

	q.X = (lhs.X*qFactor + rhs.X*otherFactor) * inverseFactor
	q.Y = (lhs.Y*qFactor + rhs.Y*otherFactor) * inverseFactor
	q.Z = (lhs.Z*qFactor + rhs.Z*otherFactor) * inverseFactor
	q.W = (lhs.W*qFactor + rhs.W*otherFactor) * inverseFactor

	return q
}

// SetSlerp performs a spherical linear interpolation at constant speed between two quaternions
// and overwrites the current quaternion with the result.
//
// lhs - The origin quaternion in the interpolation
//
// rhs - The target quaternion in the interpolation
//
// delta - The mix factor between the two quaternions- this can be any value: between 0 and 1,
// the interpolation will move linearly from the lhs to the rhs quaternion. Beyond those bounds,
// the interpolation oscillates between the two values
func (q *Quaternion[T]) SetSlerp(lhs, rhs *Quaternion[T], delta T) *Quaternion[T] {

	var rhs2 Quaternion[T]
	rhs2.SetQuaternion(rhs)

	cos := lhs.DotProduct(rhs)

	// If cos < 0, the interpolation will take the long way around the sphere
	if cos < 0 {
		rhs2.X = -rhs2.X
		rhs2.Y = -rhs2.Y
		rhs2.Z = -rhs2.Z
		rhs2.W = -rhs2.W
		cos = -cos
	}

	// Perform linear interpolation when costheta is close to 1
	if cos > 0.9999 {
		q.X = mix[T](lhs.X, rhs.X, delta)
		q.Y = mix[T](lhs.Y, rhs.Y, delta)
		q.Z = mix[T](lhs.Z, rhs.Z, delta)
		q.W = mix[T](lhs.W, rhs.W, delta)
		return q
	}

	angle := T(math.Acos(float64(cos)))
	qFactor := T(math.Sin(float64((1.0 - delta) * angle)))
	otherFactor := T(math.Sin(float64(delta * angle)))
	inverseFactor := 1.0 / T(math.Sin(float64(angle)))

	q.X = (lhs.X*qFactor + rhs.X*otherFactor) * inverseFactor
	q.Y = (lhs.Y*qFactor + rhs.Y*otherFactor) * inverseFactor
	q.Z = (lhs.Z*qFactor + rhs.Z*otherFactor) * inverseFactor
	q.W = (lhs.W*qFactor + rhs.W*otherFactor) * inverseFactor

	return q
}

// SetLerp performs a linear interpolation between two quaternions at non-constant speed and overwrites
// this quaternion with the result.
//
// lhs - The origin quaternion in the interpolation
//
// rhs - The target quaternion in the interpolation
//
// delta - The mix factor between the two quaternions and must be a value between 0 and 1
func (q *Quaternion[T]) SetLerp(lhs, rhs *Quaternion[T], delta T) *Quaternion[T] {
	q.X = lhs.X*(1-delta) + (rhs.X * delta)
	q.Y = lhs.Y*(1-delta) + (rhs.Y * delta)
	q.Z = lhs.Z*(1-delta) + (rhs.Z * delta)
	q.W = lhs.W*(1-delta) + (rhs.W * delta)

	return q.Normalize()
}

// SetMultQuaternion multiplies two quaternions together and overwrites the current contents
// of this quaternion with the results.
//
// lhs - The left operand of the multiplication operation
//
// rhs - The right operand of the multiplication operation
func (q *Quaternion[T]) SetMultQuaternion(lhs, rhs *Quaternion[T]) *Quaternion[T] {
	q.X = lhs.W*rhs.X + lhs.X*rhs.W + lhs.Y*rhs.Z - lhs.Z*rhs.Y
	q.Y = lhs.W*rhs.Y + lhs.Y*rhs.W + lhs.Z*rhs.X - lhs.X*rhs.Z
	q.Z = lhs.W*rhs.Z + lhs.Z*rhs.W + lhs.X*rhs.Y - lhs.Y*rhs.X
	q.W = lhs.W*rhs.W - lhs.X*rhs.X - lhs.Y*rhs.Y - lhs.Z*rhs.Z

	return q
}

// SetSquad interpolates between two keyframes of a SQUAD (Spherical Quadrangle Interpolation)
// path/animation. It accepts two keyframes, and the control point quaternions that correspond to those
// keyframes, and produces a spline interpolation between the two keyframes.
//
// keyframe1 - The first keyframe rotation to interpolate between
//
// keyframe2 - The second keyframe rotation to interpolate between
//
// control1 - The control point corresponding to the first keyframe. Can be produced with SetIntermediate
//
// control2 - The control point corresponding to the second keyframe. Can be produced with SetIntermediate
//
// delta - A value between 0 and 1 indicating how far to interpolate between the two keyframes
func (q *Quaternion[T]) SetSquad(keyframe1, keyframe2, control1, control2 *Quaternion[T], delta T) *Quaternion[T] {
	var qComponent, sComponent Quaternion[T]

	qComponent.SetMix(keyframe1, keyframe2, delta)
	sComponent.SetMix(control1, control2, delta)

	finalDelta := 2.0 * (1 - delta) * delta

	return q.SetMix(&qComponent, &sComponent, finalDelta)
}

// SetRotationEulers overwrites the current contents of this quaternion with quaternion
// that rotates first by yaw (y rotation), then by pitch (x rotation), and then by roll (z rotation)
//
// yawRad - Angle to rotate yaw in radians
//
// pitchRad - Angle to rotate pitch in radians
//
// rollRad - Angle to rotate roll in radians
func (q *Quaternion[T]) SetRotationEulers(rollRad, yawRad, pitchRad float64) *Quaternion[T] {
	yawCos := T(math.Cos(yawRad * 0.5))
	yawSin := T(math.Sin(yawRad * 0.5))
	pitchCos := T(math.Cos(pitchRad * 0.5))
	pitchSin := T(math.Sin(pitchRad * 0.5))
	rollCos := T(math.Cos(rollRad * 0.5))
	rollSin := T(math.Sin(rollRad * 0.5))

	q.W = pitchCos*yawCos*rollCos + pitchSin*yawSin*rollSin
	q.X = pitchSin*yawCos*rollCos - pitchCos*yawSin*rollSin
	q.Y = pitchCos*yawSin*rollCos + pitchSin*yawCos*rollSin
	q.Z = pitchCos*yawCos*rollSin - pitchSin*yawSin*rollCos

	return q
}

// SetConjugate overwrites the contents of this quaternion with the conjugate of a provided
// quaternion.
//
// other - The quaternion to conjugate
func (q *Quaternion[T]) SetConjugate(other *Quaternion[T]) *Quaternion[T] {
	q.X = -other.X
	q.Y = -other.Y
	q.Z = -other.Z
	q.W = other.W

	return q
}

// Angle retrieves the scalar angle of rotation, in radians, of the quaternion's axis-angle formulation
func (q *Quaternion[T]) Angle() T {
	if abs[T](q.W) > cosOneOverTwo[T]() {
		quatLen := T(math.Sqrt(float64(q.X*q.X + q.Y*q.Y + q.Z*q.Z)))
		return T(math.Asin(float64(quatLen))) * 2
	}

	return T(math.Acos(float64(q.W))) * 2
}

// GetAxis retrieves the vector axis of rotation of the quaternion's axis-angle formulation
//
// outAxis - A pointer to a 3-element vector to populate with the axis data
func (q *Quaternion[T]) GetAxis(outAxis *Vec3[T]) {
	tmp1 := T(1.0) - q.W*q.W

	if tmp1 <= T(0.0) {
		outAxis.X = 0
		outAxis.Y = 0
		outAxis.Z = 1

		return
	}

	tmp2 := T(1.0) / T(math.Sqrt(float64(tmp1)))

	outAxis.X = q.X * tmp2
	outAxis.Y = q.Y * tmp2
	outAxis.Z = q.Z * tmp2
}

// EulerAngles retrieves the scalar rotations, in radians, of the quaternion's
// euler angles
func (q *Quaternion[T]) EulerAngles() (yaw, pitch, roll T) {
	return q.Yaw(), q.Pitch(), q.Roll()
}

// Yaw retrieves the scalar rotation, in radians, of the quaternion's euler yaw
func (q *Quaternion[T]) Yaw() T {
	unclamped := (q.X*q.Z - q.W*q.Y) * -2
	return T(math.Asin(float64(clamp[T](unclamped, -1, 1))))
}

// Pitch retrieves the scalar rotation, in radians, of the quaternion's euler pitch
func (q *Quaternion[T]) Pitch() T {
	y := (q.Y*q.Z + q.W*q.X) * 2
	x := q.W*q.W - q.X*q.X - q.Y*q.Y + q.Z*q.Z

	if abs[T](x) < 0.0001 && abs[T](y) < 0.0001 {
		return T(math.Atan2(float64(q.X), float64(q.W))) * 2
	}

	return T(math.Atan2(float64(y), float64(x)))
}

// Roll retrieves the scalar rotation, in radians, of the quaternion's euler roll
func (q *Quaternion[T]) Roll() T {
	y := 2 * (q.X*q.Y + q.W*q.Z)
	x := q.W*q.W + q.X*q.X - q.Y*q.Y - q.Z*q.Z

	return T(math.Atan2(float64(y), float64(x)))
}

// Len retrieves the square root of the quaternion's dot product with itself
func (q *Quaternion[T]) Len() T {
	sqr := float64(q.X*q.X + q.Y*q.Y + q.Z*q.Z + q.W*q.W)
	return T(math.Sqrt(sqr))
}

// LenSqr retrieves the quaternion's dot product with itself
func (q *Quaternion[T]) LenSqr() T {
	return q.X*q.X + q.Y*q.Y + q.Z*q.Z + q.W*q.W
}

// Conjugate calculates the conjugate of the current quaternion and updates the quaternion
// to that conjugate
func (q *Quaternion[T]) Conjugate() *Quaternion[T] {
	q.X = -q.X
	q.Y = -q.Y
	q.Z = -q.Z

	return q
}

// Exp calculates the exponent of the current quaternion and updates the quaternion
// to that exponent. The exponent is the inverse of the logarithm.
func (q *Quaternion[T]) Exp() *Quaternion[T] {
	angle := T(math.Sqrt(float64(q.X*q.X + q.Y*q.Y + q.Z*q.Z)))
	if angle < 0.0001 {
		q.X = 0
		q.Y = 0
		q.Z = 0
		q.W = 1
		return q
	}

	cos := T(math.Cos(float64(angle)))
	sin := T(math.Sin(float64(angle)))
	factor := sin / angle

	q.X *= factor
	q.Y *= factor
	q.Z *= factor
	q.W = cos

	return q
}

// Log calculates the logarithm of the current quaternion and updates the quaternion
// to that logarithm. The logarithm is the inverse of the exponent.
func (q *Quaternion[T]) Log() *Quaternion[T] {
	lenSquared := q.LenSqr()
	quatLen := T(math.Sqrt(float64(lenSquared)))

	if lenSquared < 0.0001 {
		if q.W > 0 {
			q.X = 0
			q.Y = 0
			q.Z = 0
			q.W = T(math.Log(float64(q.W)))
			return q
		} else if q.W < 0 {
			q.X = math.Pi
			q.Y = 0
			q.Z = 0
			q.W = T(math.Log(float64(-q.W)))
			return q
		}

		q.X = T(math.Inf(1))
		q.Y = T(math.Inf(1))
		q.Z = T(math.Inf(1))
		q.W = T(math.Inf(1))
		return q
	}

	t := T(math.Atan2(float64(quatLen), float64(q.W))) / quatLen
	quatLen2 := lenSquared + q.W*q.W

	q.X *= t
	q.Y *= t
	q.Z *= t
	q.W = T(math.Log(float64(quatLen2))) * 0.5

	return q
}

// Pow raises the current quaternion to a provided power and updates the quaternion to the
// calculated quaternion. Pow is understood in terms of the quaternion multiplication operation,
// so squaring a quaternion is the same as multiplying a quaternion against itself, cubing a
// quaternion is the same as multiplying it against its square, and so forth.
//
// y - A scalar value indicating the power to raise the quaternion to. Can be any value, a quaternion
// raised to the 0th power is the identity quaternion.
func (q *Quaternion[T]) Pow(y T) *Quaternion[T] {
	if abs[T](y) < 0.0001 {
		q.X = 0
		q.Y = 0
		q.Z = 0
		q.W = 1
		return q
	}

	magnitude := q.Len()

	var angle T
	if abs[T](q.W/magnitude) > cosOneOverTwo[T]() {
		vectorMagnitude := q.X*q.X + q.Y*q.Y + q.Z*q.Z

		if vectorMagnitude < 0.0001 {
			q.X = 0
			q.Y = 0
			q.Z = 0
			q.W = T(math.Pow(float64(q.W), float64(y)))
			return q
		}

		angle = T(math.Asin(math.Sqrt(float64(vectorMagnitude)) / float64(magnitude)))
	} else {
		angle = T(math.Acos(float64(q.W / magnitude)))
	}

	newAngle := angle * y
	div := T(math.Sin(float64(newAngle)) / math.Sin(float64(angle)))
	mag := T(math.Pow(float64(magnitude), float64(y-1)))

	q.X = q.X * div * mag
	q.Y = q.Y * div * mag
	q.Z = q.Z * div * mag
	q.W = T(math.Cos(float64(newAngle))) * magnitude * mag

	return q
}

// Normalize updates the current quaternion to be a unit quaternion
func (q *Quaternion[T]) Normalize() *Quaternion[T] {
	oneOverLen := 1.0 / q.Len()

	q.X *= oneOverLen
	q.Y *= oneOverLen
	q.Z *= oneOverLen
	q.W *= oneOverLen

	return q
}

// Inverse updates the current quaternion to be its own inverse. For unit quaternions,
// this is the same as the conjugate, but is less performant in order to cover non-unit cases.
func (q *Quaternion[T]) Inverse() *Quaternion[T] {
	inverseDotProduct := 1.0 / (q.X*q.X + q.Y*q.Y + q.Z*q.Z + q.W*q.W)
	q.X = -q.X * inverseDotProduct
	q.Y = -q.Y * inverseDotProduct
	q.Z = -q.Z * inverseDotProduct
	q.W = q.W * inverseDotProduct

	return q
}

// Mix performs an oriented spherical linear interpolation at constant speed between this quaternion
// and another provided quaternion and updates this quaternion to the result
//
// other - The target quaternion in the interpolation
//
// delta - The mix factor between the two quaternions- this can be any value: between 0 and 1,
// the interpolation will move linearly from this quaternion to the other quaternion. Beyond those bounds,
// the interpolation oscillates between the two values
func (q *Quaternion[T]) Mix(other *Quaternion[T], delta T) *Quaternion[T] {
	cosTheta := q.DotProduct(other)

	// Perform a linear interpolation when cosTheta is close to 1 to avoid side effect of sin(angle) becoming a zero denominator
	if cosTheta > 0.9999 {
		q.X = mix[T](q.X, other.X, delta)
		q.Y = mix[T](q.Y, other.Y, delta)
		q.Z = mix[T](q.Z, other.Z, delta)
		q.W = mix[T](q.W, other.W, delta)

		return q
	}

	angle := math.Acos(float64(cosTheta))

	//(sin((static_cast<T>(1) - a) * angle) * x + sin(a * angle) * y) / sin(angle);
	qFactor := T(math.Sin(float64(1.0-delta) * angle))
	otherFactor := T(math.Sin(float64(delta) * angle))
	inverseFactor := 1.0 / T(math.Sin(float64(angle)))

	q.X = (q.X*qFactor + other.X*otherFactor) * inverseFactor
	q.Y = (q.Y*qFactor + other.Y*otherFactor) * inverseFactor
	q.Z = (q.Z*qFactor + other.Z*otherFactor) * inverseFactor
	q.W = (q.W*qFactor + other.W*otherFactor) * inverseFactor

	return q
}

// Slerp performs a spherical linear interpolation at constant speed between this quaternion and another
// provided quaternion and updates this quaternion with the result
//
// other - The target quaternion in the interpolation
//
// delta - The mix factor between the two quaternions- this can be any value: between 0 and 1,
// the interpolation will move linearly from this quaternion to the other quaternion. Beyond those bounds,
// the interpolation oscillates between the two values
func (q *Quaternion[T]) Slerp(other *Quaternion[T], delta T) *Quaternion[T] {
	cosTheta := q.DotProduct(other)

	var tmpOther Quaternion[T]
	tmpOther.SetQuaternion(other)

	// If cosTheta < 0, the interpolation will take the long way around the sphere
	if cosTheta < 0 {
		tmpOther.X = -tmpOther.X
		tmpOther.Y = -tmpOther.Y
		tmpOther.Z = -tmpOther.Z
		tmpOther.W = -tmpOther.W
		cosTheta = -cosTheta
	}

	// Perform linear interpolation when costheta is close to 1
	if cosTheta > 0.9999 {
		q.X = mix[T](q.X, other.X, delta)
		q.Y = mix[T](q.Y, other.Y, delta)
		q.Z = mix[T](q.Z, other.Z, delta)
		q.W = mix[T](q.W, other.W, delta)
		return q
	}

	angle := T(math.Acos(float64(cosTheta)))
	qFactor := T(math.Sin(float64((1.0 - delta) * angle)))
	otherFactor := T(math.Sin(float64(delta * angle)))
	inverseFactor := 1.0 / T(math.Sin(float64(angle)))

	q.X = (q.X*qFactor + other.X*otherFactor) * inverseFactor
	q.Y = (q.Y*qFactor + other.Y*otherFactor) * inverseFactor
	q.Z = (q.Z*qFactor + other.Z*otherFactor) * inverseFactor
	q.W = (q.W*qFactor + other.W*otherFactor) * inverseFactor

	return q
}

// Lerp performs a linear interpolation between this quaternion and another provided quaternion
// at non-constant speed and updates this quaternion with the result.
//
// other - The target quaternion in the interpolation
//
// delta - The mix factor between the two quaternions, must be a value between 0 and 1
func (q *Quaternion[T]) Lerp(other *Quaternion[T], delta T) *Quaternion[T] {
	q.X = q.X*(1-delta) + (other.X * delta)
	q.Y = q.Y*(1-delta) + (other.Y * delta)
	q.Z = q.Z*(1-delta) + (other.Z * delta)
	q.W = q.W*(1-delta) + (other.W * delta)

	return q.Normalize()
}

// RotateX applies a rotation to this quaternion that rotates around the x axis by the specified amount
//
// angleRad - The angle to rotate around the x axis in radians
func (q *Quaternion[T]) RotateX(angleRad float64) *Quaternion[T] {
	sin := T(math.Sin(angleRad * 0.5))
	cos := T(math.Cos(angleRad * 0.5))

	x := q.W*sin + q.X*cos
	y := q.Y*cos + q.Z*sin
	z := q.Z*cos - q.Y*sin
	w := q.W*cos - q.X*sin

	q.X = x
	q.Y = y
	q.Z = z
	q.W = w

	return q
}

// RotateY applies a rotation to this quaternion that rotates around the y axis by the specified amount
//
// angleRad - The angle to rotate around the y axis in radians
func (q *Quaternion[T]) RotateY(angleRad float64) *Quaternion[T] {
	sin := T(math.Sin(angleRad * 0.5))
	cos := T(math.Cos(angleRad * 0.5))

	x := q.X*cos - q.Z*sin
	y := q.W*sin + q.Y*cos
	z := q.Z*cos + q.X*sin
	w := q.W*cos - q.Y*sin

	q.X = x
	q.Y = y
	q.Z = z
	q.W = w

	return q
}

// RotateZ applies a rotation to this quaternion that rotates around the z axis by the specified amount
//
// angleRad - The angle to rotate around the z axis in radians
func (q *Quaternion[T]) RotateZ(angleRad float64) *Quaternion[T] {
	sin := T(math.Sin(angleRad * 0.5))
	cos := T(math.Cos(angleRad * 0.5))

	x := q.X*cos + q.Y*sin
	y := q.Y*cos - q.X*sin
	z := q.W*sin + q.Z*cos
	w := q.W*cos - q.Z*sin

	q.X = x
	q.Y = y
	q.Z = z
	q.W = w

	return q
}

// RotateAroundAxis rotates this quaternion by applying a rotation around the provided axis by the
// provided angle in radians
//
// axis - A 3-element vector that is normal to the angle of rotation. It does not need to be normalized.
//
// angleRad - The amount to rotate in radians
func (q *Quaternion[T]) RotateAroundAxis(axis *Vec3[T], angleRad float64) *Quaternion[T] {
	var unitAxis Vec3[T]
	unitAxis.SetVec3(axis).Normalize()

	sin := T(math.Sin(angleRad * 0.5))

	angleAxisQuat := Quaternion[T]{
		X: unitAxis.X * sin,
		Y: unitAxis.Y * sin,
		Z: unitAxis.Z * sin,
		W: T(math.Cos(angleRad * 0.5)),
	}

	return q.MultQuaternion(&angleAxisQuat)
}

// MultQuaternion multiplies this quaternion with a provided quaternion and updates this quaternion
// with the result
//
// other - The quaternion to use as the right operand in the multiplication operation
func (q *Quaternion[T]) MultQuaternion(other *Quaternion[T]) *Quaternion[T] {
	x := q.W*other.X + q.X*other.W + q.Y*other.Z - q.Z*other.Y
	y := q.W*other.Y + q.Y*other.W + q.Z*other.X - q.X*other.Z
	z := q.W*other.Z + q.Z*other.W + q.X*other.Y - q.Y*other.X
	w := q.W*other.W - q.X*other.X - q.Y*other.Y - q.Z*other.Z

	q.X = x
	q.Y = y
	q.Z = z
	q.W = w

	return q
}

// DotProduct calculates and returns the dot product of this quaternion with another provided quaternion
//
// other - The quaternion to use as the right operand in the dot product operation
func (q *Quaternion[T]) DotProduct(other *Quaternion[T]) T {
	return q.X*other.X + q.Y*other.Y + q.Z*other.Z + q.W*other.W
}

// Equal returns true if every entry in this quaternion is equal to every entry in the provided quaternion
//
// other - The quaternion to compare this quaternion to
//
// epsilon - The epsilon value to use in floating point comparisons. This much floating point
// drift is permitted before the method returns false. 0.0001 is a common epsilon value
func (q *Quaternion[T]) Equal(other *Quaternion[T], epsilon T) bool {
	xDiff := abs[T](q.X - other.X)
	yDiff := abs[T](q.Y - other.Y)
	zDiff := abs[T](q.Z - other.Z)
	wDiff := abs[T](q.W - other.W)

	return xDiff < epsilon && yDiff < epsilon && zDiff < epsilon && wDiff < epsilon
}
