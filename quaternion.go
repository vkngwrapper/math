package math

import "math"

type Quaternion[T FloatingPoint] struct {
	X, Y, Z, W T
}

func (q *Quaternion[T]) SetIdentity() *Quaternion[T] {
	q.X = 0
	q.Y = 0
	q.Z = 0
	q.W = 1

	return q
}

func (q *Quaternion[T]) SetAxisAngle(axis *Vec3[T], angle T) *Quaternion[T] {
	s := T(math.Sin(float64(angle * 0.5)))
	q.X = axis.X * s
	q.Y = axis.Y * s
	q.Z = axis.Z * s
	q.W = angle * 0.5

	return q
}

func (q *Quaternion[T]) SetVectorRotation(origin *Vec3[T], destination *Vec3[T]) *Quaternion[T] {
	cosTheta := origin.DotProduct(destination)

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
		q.W = T(math.Pi)
		return q
	}

	var rotationAxis Vec3[T]
	rotationAxis.SetCrossProduct(origin, destination)

	lengthSquared := (1.0 + cosTheta) * 2.0
	length := T(math.Sqrt(float64(lengthSquared)))
	inverse := 1.0 / length

	q.X = rotationAxis.X * inverse
	q.Y = rotationAxis.Y * inverse
	q.Z = rotationAxis.Z * inverse
	q.W = length * 0.5

	return q
}

func (q *Quaternion[T]) SetQuaternion(other *Quaternion[T]) *Quaternion[T] {
	q.X = other.X
	q.Y = other.Y
	q.Z = other.Z
	q.W = other.W

	return q
}

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

func (q *Quaternion[T]) SetIntermediate(prev *Quaternion[T], current *Quaternion[T], next *Quaternion[T]) *Quaternion[T] {
	var inverseCurrent Quaternion[T]
	inverseCurrent.Inverse()

	var nextDiff Quaternion[T]
	nextDiff.SetMultQuaternion(next, &inverseCurrent).Log()

	var prevDiff Quaternion[T]
	prevDiff.SetMultQuaternion(prev, &inverseCurrent).Log()

	factor := T(-1.0 / 4.0)
	preExp := Quaternion[T]{
		X: (nextDiff.X + prevDiff.X) * factor,
		Y: (nextDiff.Y + prevDiff.Y) * factor,
		Z: (nextDiff.Z + prevDiff.Z) * factor,
		W: (nextDiff.W + prevDiff.W) * factor,
	}
	return preExp.Exp().MultQuaternion(current)
}

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

	angle := T(math.Acos(float64(cosTheta)))
	qFactor := T(math.Sin(float64(1.0-delta))) * angle
	otherFactor := T(math.Sin(float64(delta * angle)))
	inverseFactor := 1.0 / T(math.Sin(float64(angle)))

	q.X = (lhs.X*qFactor + rhs.X*otherFactor) * inverseFactor
	q.Y = (lhs.Y*qFactor + rhs.Y*otherFactor) * inverseFactor
	q.Z = (lhs.Z*qFactor + rhs.Z*otherFactor) * inverseFactor
	q.W = (lhs.W*qFactor + rhs.W*otherFactor) * inverseFactor

	return q
}

func (q *Quaternion[T]) SetShortMix(lhs, rhs *Quaternion[T], delta T) *Quaternion[T] {
	if delta <= 0 {
		return q
	} else if delta >= 1 {
		q.X = rhs.X
		q.Y = rhs.Y
		q.Z = rhs.Z
		q.W = rhs.W
		return q
	}

	var rhs2 Quaternion[T]
	rhs2.SetQuaternion(rhs)

	cos := lhs.DotProduct(rhs)

	if cos < 0 {
		rhs2.X = -rhs2.X
		rhs2.Y = -rhs2.Y
		rhs2.Z = -rhs2.Z
		rhs2.W = -rhs2.W
		cos = -cos
	}

	var k0, k1 T
	if cos > 0.9999 {
		k0 = 1.0 - delta
		k1 = delta
	} else {
		sin := T(math.Sqrt(float64(1.0 - cos*cos)))
		angle := T(math.Atan2(float64(sin), float64(cos)))
		oneOverSin := 1.0 / sin
		k0 = T(math.Sin(float64((1.0-delta)*angle))) * oneOverSin
		k1 = T(math.Sin(float64(delta*angle))) * oneOverSin
	}

	q.X = k0*lhs.X + k1*rhs2.X
	q.Y = k0*lhs.Y + k1*rhs2.Y
	q.Z = k0*lhs.Z + k1*rhs2.Z
	q.W = k0*lhs.W + k1*rhs2.W

	return q
}

func (q *Quaternion[T]) SetSlerp(lhs, rhs *Quaternion[T], delta T) *Quaternion[T] {
	cosTheta := lhs.DotProduct(rhs)

	var tmpRHS Quaternion[T]
	tmpRHS.SetQuaternion(rhs)

	// If cosTheta < 0, the interpolation will take the long way around the sphere
	if cosTheta < 0 {
		tmpRHS.X = -tmpRHS.X
		tmpRHS.Y = -tmpRHS.Y
		tmpRHS.Z = -tmpRHS.Z
		tmpRHS.W = -tmpRHS.W
		cosTheta = -cosTheta
	}

	// Perform linear interpolation when costheta is close to 1
	if cosTheta > 0.9999 {
		q.X = mix[T](lhs.X, rhs.X, delta)
		q.Y = mix[T](lhs.Y, rhs.Y, delta)
		q.Z = mix[T](lhs.Z, rhs.Z, delta)
		q.W = mix[T](lhs.W, rhs.W, delta)
		return q
	}

	angle := T(math.Acos(float64(cosTheta)))
	qFactor := T(math.Sin(float64((1.0 - delta) * angle)))
	otherFactor := T(math.Sin(float64(delta * angle)))
	inverseFactor := 1.0 / T(math.Sin(float64(angle)))

	q.X = (lhs.X*qFactor + rhs.X*otherFactor) * inverseFactor
	q.Y = (lhs.Y*qFactor + rhs.Y*otherFactor) * inverseFactor
	q.Z = (lhs.Z*qFactor + rhs.Z*otherFactor) * inverseFactor
	q.W = (lhs.W*qFactor + rhs.W*otherFactor) * inverseFactor

	return q
}

func (q *Quaternion[T]) SetLerp(lhs, rhs *Quaternion[T], delta T) *Quaternion[T] {
	q.X = lhs.X*(1-delta) + (rhs.X * delta)
	q.Y = lhs.Y*(1-delta) + (rhs.Y * delta)
	q.Z = lhs.Z*(1-delta) + (rhs.Z * delta)
	q.W = lhs.W*(1-delta) + (rhs.W * delta)

	return q.Normalize()
}

func (q *Quaternion[T]) SetMultQuaternion(lhs, rhs *Quaternion[T]) *Quaternion[T] {
	x := lhs.W*rhs.X + lhs.X*rhs.W + lhs.Y*rhs.Z - lhs.Z*rhs.Y
	y := lhs.W*rhs.Y + lhs.Y*rhs.W + lhs.Z*rhs.X - lhs.X*rhs.Z
	z := lhs.W*rhs.Z + lhs.Z*rhs.W + lhs.X*rhs.Y - lhs.Y*rhs.X
	w := lhs.W*rhs.W - lhs.X*rhs.X - lhs.Y*rhs.Y - lhs.Z*rhs.Z

	q.X = x
	q.Y = y
	q.Z = z
	q.W = w

	return q
}

func (q *Quaternion[T]) SetQuad(q1, q2, s1, s2 *Quaternion[T], delta T) *Quaternion[T] {
	var qComponent, sComponent Quaternion[T]

	qComponent.SetMix(q1, q2, delta)
	sComponent.SetMix(s1, s2, delta)

	finalDelta := 2.0 * (1 - delta) * delta

	return q.SetMix(&qComponent, &sComponent, finalDelta)
}

func (q *Quaternion[T]) Angle() T {
	if abs[T](q.W) > cosOneOverTwo[T]() {
		quatLen := T(math.Sqrt(float64(q.X*q.X + q.Y*q.Y + q.Z*q.Z)))
		return T(math.Asin(float64(quatLen))) * 2
	}

	return T(math.Acos(float64(q.W))) * 2
}

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

func (q *Quaternion[T]) EulerAngles() (yaw, pitch, roll T) {
	return q.Yaw(), q.Pitch(), q.Roll()
}

func (q *Quaternion[T]) Yaw() T {
	unclamped := (q.X*q.Z - q.W*q.Y) * -2
	return T(math.Asin(float64(clamp[T](unclamped, -1, 1))))
}

func (q *Quaternion[T]) Pitch() T {
	y := (q.Y*q.Z + q.W*q.X) * 2
	x := q.W*q.W - q.X*q.X - q.Y*q.Y + q.Z*q.Z

	if abs[T](x) < 0.0001 && abs[T](y) < 0.0001 {
		return T(math.Atan2(float64(q.X), float64(q.W))) * 2
	}

	return T(math.Atan2(float64(y), float64(x)))
}

func (q *Quaternion[T]) Roll() T {
	y := 2 * (q.X*q.Y + q.W*q.Z)
	x := q.W*q.W + q.X*q.X - q.Y*q.Y - q.Z*q.Z

	return T(math.Atan2(float64(y), float64(x)))
}

func (q *Quaternion[T]) ExtractRealComponent() T {
	w := 1.0 - q.X*q.X - q.Y*q.Y - q.Z*q.Z

	if w < 0 {
		return 0
	} else {
		return -T(math.Sqrt(float64(w)))
	}
}

func (q *Quaternion[T]) Len() T {
	sqr := float64(q.X*q.X + q.Y*q.Y + q.Z*q.Z + q.W*q.W)
	return T(math.Sqrt(sqr))
}

func (q *Quaternion[T]) LenSqr() T {
	return q.X*q.X + q.Y*q.Y + q.Z*q.Z + q.W*q.W
}

func (q *Quaternion[T]) Conjugate() *Quaternion[T] {
	q.X = -q.X
	q.Y = -q.Y
	q.Z = -q.Z

	return q
}

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

func (q *Quaternion[T]) Normalize() *Quaternion[T] {
	oneOverLen := 1.0 / q.Len()

	q.X *= oneOverLen
	q.Y *= oneOverLen
	q.Z *= oneOverLen
	q.W *= oneOverLen

	return q
}

func (q *Quaternion[T]) Inverse() *Quaternion[T] {
	inverseDotProduct := 1.0 / (q.X*q.X + q.Y*q.Y + q.Z*q.Z + q.W*q.W)
	q.X = -q.X * inverseDotProduct
	q.Y = -q.Y * inverseDotProduct
	q.Z = -q.Z * inverseDotProduct
	q.W = q.W * inverseDotProduct

	return q
}

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

	angle := T(math.Acos(float64(cosTheta)))
	qFactor := T(math.Sin(float64(1.0-delta))) * angle
	otherFactor := T(math.Sin(float64(delta * angle)))
	inverseFactor := 1.0 / T(math.Sin(float64(angle)))

	q.X = (q.X*qFactor + other.X*otherFactor) * inverseFactor
	q.Y = (q.Y*qFactor + other.Y*otherFactor) * inverseFactor
	q.Z = (q.Z*qFactor + other.Z*otherFactor) * inverseFactor
	q.W = (q.W*qFactor + other.W*otherFactor) * inverseFactor

	return q
}

func (q *Quaternion[T]) ShortMix(other *Quaternion[T], delta T) *Quaternion[T] {
	if delta <= 0 {
		return q
	} else if delta >= 1 {
		q.X = other.X
		q.Y = other.Y
		q.Z = other.Z
		q.W = other.W
		return q
	}

	var other2 Quaternion[T]
	other2.SetQuaternion(other)

	cos := q.DotProduct(other)

	if cos < 0 {
		other2.X = -other2.X
		other2.Y = -other2.Y
		other2.Z = -other2.Z
		other2.W = -other2.W
		cos = -cos
	}

	var k0, k1 T
	if cos > 0.9999 {
		k0 = 1.0 - delta
		k1 = delta
	} else {
		sin := T(math.Sqrt(float64(1.0 - cos*cos)))
		angle := T(math.Atan2(float64(sin), float64(cos)))
		oneOverSin := 1.0 / sin
		k0 = T(math.Sin(float64((1.0-delta)*angle))) * oneOverSin
		k1 = T(math.Sin(float64(delta*angle))) * oneOverSin
	}

	q.X = k0*q.X + k1*other2.X
	q.Y = k0*q.Y + k1*other2.Y
	q.Z = k0*q.Z + k1*other2.Z
	q.W = k0*q.W + k1*other2.W

	return q
}

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

func (q *Quaternion[T]) Lerp(other *Quaternion[T], delta T) *Quaternion[T] {
	q.X = q.X*(1-delta) + (other.X * delta)
	q.Y = q.Y*(1-delta) + (other.Y * delta)
	q.Z = q.Z*(1-delta) + (other.Z * delta)
	q.W = q.W*(1-delta) + (other.W * delta)

	return q.Normalize()
}

func (q *Quaternion[T]) RotateNormalizedAxis(unitAxis *Vec3[T], angle T) *Quaternion[T] {
	sin := T(math.Sin(float64(angle * 0.5)))

	angleAxisQuat := Quaternion[T]{
		X: unitAxis.X * sin,
		Y: unitAxis.Y * sin,
		Z: unitAxis.Z * sin,
		W: T(math.Cos(float64(angle * 0.5))),
	}

	return q.MultQuaternion(&angleAxisQuat)
}

func (q *Quaternion[T]) Rotate(axis *Vec3[T], angle T) *Quaternion[T] {
	axisLen := axis.Len()

	var unitAxis Vec3[T]
	unitAxis.SetVec3(axis)

	if abs[T](axisLen-1.0) > 0.0001 {
		oneOverLen := 1.0 / axisLen
		unitAxis.X *= oneOverLen
		unitAxis.Y *= oneOverLen
		unitAxis.Z *= oneOverLen
	}

	return q.RotateNormalizedAxis(&unitAxis, angle)
}

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

func (q *Quaternion[T]) DotProduct(other *Quaternion[T]) T {
	return q.X*other.X + q.Y*other.Y + q.Z*other.Z + q.W*other.W
}

func (q *Quaternion[T]) Equal(other *Quaternion[T], epsilon T) bool {
	xDiff := abs[T](q.X - other.X)
	yDiff := abs[T](q.Y - other.Y)
	zDiff := abs[T](q.Z - other.Z)
	wDiff := abs[T](q.W - other.W)

	return xDiff < epsilon && yDiff < epsilon && zDiff < epsilon && wDiff < epsilon
}

func (q *Quaternion[T]) GreaterThan(other *Quaternion[T]) bool {
	return q.X > other.X && q.Y > other.Y && q.Z > other.Z && q.W > other.W
}

func (q *Quaternion[T]) GreaterThanEqual(other *Quaternion[T]) bool {
	return q.X >= other.X && q.Y >= other.Y && q.Z >= other.Z && q.W >= other.W
}

func (q *Quaternion[T]) LessThan(other *Quaternion[T]) bool {
	return q.X < other.X && q.Y < other.Y && q.Z < other.Z && q.W < other.W
}

func (q *Quaternion[T]) LessThanEqual(other *Quaternion[T]) bool {
	return q.X <= other.X && q.Y <= other.Y && q.Z <= other.Z && q.W <= other.W
}
