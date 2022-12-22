package math

import "math"

type Mat3Row[T FloatingPoint] [3]T

// Mat3x3 is a two-dimensional array of floating point-compatible values that can be used
// for 3x3 matrix arithmetic and affine 3d matrix transformations
type Mat3x3[T FloatingPoint] [3]Mat3Row[T]

// SetIdentity overwrites the current contents of the matrix with the identity matrix
func (m *Mat3x3[T]) SetIdentity() {
	m[0][0] = 1
	m[0][1] = 0
	m[0][2] = 0
	m[1][0] = 0
	m[1][1] = 1
	m[1][2] = 0
	m[2][0] = 0
	m[2][1] = 0
	m[2][2] = 1
}

// SetDiagonalScalar overwrites the current contents of the matrix with an identity matrix,
// multiplied by the provided scalar value
//
// s - The scalar value to multiply the matrix entries by
func (m *Mat3x3[T]) SetDiagonalScalar(s T) {
	m[0][0] = s
	m[0][1] = 0
	m[0][2] = 0
	m[1][0] = 0
	m[1][1] = s
	m[1][2] = 0
	m[2][0] = 0
	m[2][1] = 0
	m[2][2] = s
}

// SetDiagonalVector overwrites the current contents of the matrix with an identity matrix,
// in which the 1 values are replaced with the contents of a provided vector, with each vector
// element corresponding to a column in the matrix
//
// v - The vector containing the diagonal entries that will be populated into the matrix
func (m *Mat3x3[T]) SetDiagonalVector(v *Vec3[T]) {
	m[0][0] = v.X
	m[0][1] = 0
	m[0][2] = 0
	m[1][0] = 0
	m[1][1] = v.Y
	m[1][2] = 0
	m[2][0] = 0
	m[2][1] = 0
	m[2][2] = v.Z
}

// SetQuaternion overwrites the current contents of the matrix with a rotation matrix
// providing the same rotation as a provided quaternion
//
// other - A quaternion whose rotation will be written into this matrix
func (m *Mat3x3[T]) SetQuaternion(other *Quaternion[T]) {
	xx := other.X * other.X
	yy := other.Y * other.Y
	zz := other.Z * other.Z
	xz := other.X * other.Z
	xy := other.X * other.Y
	yz := other.Y * other.Z
	wx := other.W * other.X
	wy := other.W * other.Y
	wz := other.W * other.Z

	m[0][0] = 1.0 - 2.0*(yy+zz)
	m[0][1] = 2.0 * (xy + wz)
	m[0][2] = 2.0 * (xz - wy)
	m[1][0] = 2.0 * (xy - wz)
	m[1][1] = 1.0 - 2.0*(xx+zz)
	m[1][2] = 2.0 * (yz + wx)
	m[2][0] = 2.0 * (xz + wy)
	m[2][1] = 2.0 * (yz - wx)
	m[2][2] = 1.0 - 2.0*(xx+yy)
}

// SetMat2x2 overwrites the current contents of the matrix with the identity matrix and
// layers the provided 2x2 matrix over the upper left 2x2 area of the matrix
//
// other - The 2x2 matrix to initialize from
func (m *Mat3x3[T]) SetMat2x2(other *Mat3x3[T]) {
	m[0][0] = other[0][0]
	m[0][1] = other[0][1]
	m[0][2] = 0
	m[1][0] = other[1][0]
	m[1][1] = other[1][1]
	m[1][2] = 0
	m[2][0] = 0
	m[2][1] = 0
	m[2][2] = 1
}

// SetMat3x3 overwrites the current contents of the matrix with the contents of the provided
// matrix
//
// other - The matrix to initialize from
func (m *Mat3x3[T]) SetMat3x3(other *Mat3x3[T]) {
	m[0][0] = other[0][0]
	m[0][1] = other[0][1]
	m[0][2] = other[0][2]
	m[1][0] = other[1][0]
	m[1][1] = other[1][1]
	m[1][2] = other[1][2]
	m[2][0] = other[2][0]
	m[2][1] = other[2][1]
	m[2][2] = other[2][2]
}

// SetMat4x4 overwrites the current contents of the matrix with the upper left
// 3x3 area of a 4x4 matrix
//
// other - The 4x4 matrix to initialize from
func (m *Mat3x3[T]) SetMat4x4(other *Mat4x4[T]) {
	m[0][0] = other[0][0]
	m[0][1] = other[0][1]
	m[0][2] = other[0][2]
	m[1][0] = other[1][0]
	m[1][1] = other[1][1]
	m[1][2] = other[1][2]
	m[2][0] = other[2][0]
	m[2][1] = other[2][1]
	m[2][2] = other[2][2]
}

// SetColMajor overwrites the current contents of the matrix with 3 matrix columns
// passed in as 3-element vectors
//
// c1 - The first column of the matrix
// c2 - The second column of the matrix
// c3 - The third column of the matrix
func (m *Mat3x3[T]) SetColMajor(c1, c2, c3 *Vec3[T]) {
	m[0][0] = c1.X
	m[0][1] = c1.Y
	m[0][2] = c1.Z
	m[1][0] = c2.X
	m[1][1] = c2.Y
	m[1][2] = c2.Z
	m[2][0] = c3.X
	m[2][1] = c3.Y
	m[2][2] = c3.Z
}

// SetRowMajor overwrites the current contents of the matrix with 3 matrix rows
// passed in as 3-element vectors
//
// r1 - The first row of the matrix
// r2 - The second row of the matrix
// r3 - The third row of the matrix
func (m *Mat3x3[T]) SetRowMajor(r1, r2, r3 *Vec3[T]) {
	m[0][0] = r1.X
	m[1][0] = r1.Y
	m[2][0] = r2.Z
	m[0][1] = r2.X
	m[1][1] = r2.Y
	m[2][1] = r2.Z
	m[0][2] = r3.X
	m[1][2] = r3.Y
	m[2][2] = r3.Z
}

// SetMultMat3x3 multiplies two 3x3 matrices together and overwrites the current contents
// of this matrix with the results. In Matrix Multiplication between transform matrices,
// the left matrix is unintuitively applied to the right matrix rather than the other way around.
// You may prefer to use SetApplyTransform.
//
// lhs - The left operand of the multiplication operation
// rhs - The right operand of the multiplication operation
func (m *Mat3x3[T]) SetMultMat3x3(lhs, rhs *Mat3x3[T]) {
	m00 := lhs[0][0]*rhs[0][0] + lhs[1][0]*rhs[0][1] + lhs[2][0]*rhs[0][2]
	m10 := lhs[0][0]*rhs[1][0] + lhs[1][0]*rhs[1][1] + lhs[2][0]*rhs[1][2]
	m20 := lhs[0][0]*rhs[2][0] + lhs[1][0]*rhs[2][1] + lhs[2][0]*rhs[2][2]

	m01 := lhs[0][1]*rhs[0][0] + lhs[1][1]*rhs[0][1] + lhs[2][1]*rhs[0][2]
	m11 := lhs[0][1]*rhs[1][0] + lhs[1][1]*rhs[1][1] + lhs[2][1]*rhs[1][2]
	m21 := lhs[0][1]*rhs[2][0] + lhs[1][1]*rhs[2][1] + lhs[2][1]*rhs[2][2]

	m02 := lhs[0][2]*rhs[0][0] + lhs[1][2]*rhs[0][1] + lhs[2][2]*rhs[0][2]
	m12 := lhs[0][2]*rhs[1][0] + lhs[1][2]*rhs[1][1] + lhs[2][2]*rhs[1][2]
	m22 := lhs[0][2]*rhs[2][0] + lhs[1][2]*rhs[2][1] + lhs[2][2]*rhs[2][2]

	m[0][0] = m00
	m[0][1] = m01
	m[0][2] = m02
	m[1][0] = m10
	m[1][1] = m11
	m[1][2] = m12
	m[2][0] = m20
	m[2][1] = m21
	m[2][2] = m22
}

// SetApplyTransform applies the right transform matrix to the left transform matrix and overwrites
// the current contents of this matrix with the result.  "ApplyTransform" is just a matrix multiply
// with the operands reversed.
//
// lhs - The matrix having a transform applied to it
//
// rhs - The transform being applied
func (m *Mat3x3[T]) SetApplyTransform(lhs, rhs *Mat3x3[T]) {
	m00 := rhs[0][0]*lhs[0][0] + rhs[1][0]*lhs[0][1] + rhs[2][0]*lhs[0][2]
	m10 := rhs[0][0]*lhs[1][0] + rhs[1][0]*lhs[1][1] + rhs[2][0]*lhs[1][2]
	m20 := rhs[0][0]*lhs[2][0] + rhs[1][0]*lhs[2][1] + rhs[2][0]*lhs[2][2]

	m01 := rhs[0][1]*lhs[0][0] + rhs[1][1]*lhs[0][1] + rhs[2][1]*lhs[0][2]
	m11 := rhs[0][1]*lhs[1][0] + rhs[1][1]*lhs[1][1] + rhs[2][1]*lhs[1][2]
	m21 := rhs[0][1]*lhs[2][0] + rhs[1][1]*lhs[2][1] + rhs[2][1]*lhs[2][2]

	m02 := rhs[0][2]*lhs[0][0] + rhs[1][2]*lhs[0][1] + rhs[2][2]*lhs[0][2]
	m12 := rhs[0][2]*lhs[1][0] + rhs[1][2]*lhs[1][1] + rhs[2][2]*lhs[1][2]
	m22 := rhs[0][2]*lhs[2][0] + rhs[1][2]*lhs[2][1] + rhs[2][2]*lhs[2][2]

	m[0][0] = m00
	m[0][1] = m01
	m[0][2] = m02
	m[1][0] = m10
	m[1][1] = m11
	m[1][2] = m12
	m[2][0] = m20
	m[2][1] = m21
	m[2][2] = m22
}

// GetAxisAngle retrieves the axis and angle of rotation for this matrix, if it is a rotation matrix.
// If it is not a rotation matrix, the results are undefined.
//
// outAxis - A pointer to a 3-element vector that will be populated with a unit vector normal
// to the angle of rotation.
//
// outAngleRad - A pointer to a float64 that will be populated with the amount to rotate in radians
func (m *Mat3x3[T]) GetAxisAngle(outAxis *Vec3[T], outAngleRad *float64) {
	epsilon := T(0.01)
	epsilon2 := T(0.1)

	if (abs[T](m[1][0]-m[0][1]) < epsilon) &&
		(abs[T](m[2][0]-m[0][2]) < epsilon) &&
		(abs[T](m[2][1]-m[1][2]) < epsilon) {

		if (abs[T](m[1][0]+m[0][1]) < epsilon2) &&
			(abs[T](m[2][0]+m[0][2]) < epsilon2) &&
			(abs[T](m[2][1]+m[1][2]) < epsilon2) &&
			(abs[T](m[0][0]+m[1][1]+m[2][2]-T(3)) < epsilon2) {

			*outAngleRad = 0
			outAxis.X = T(1)
			outAxis.Y = T(0)
			outAxis.Z = T(0)
			return
		}

		*outAngleRad = math.Pi
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
	if angleCos < -1 {
		*outAngleRad = math.Pi
	} else if angleCos > 1 {
		*outAngleRad = 0
	} else {
		*outAngleRad = math.Acos(float64(angleCos))
	}

	outAxis.X = (m[1][2] - m[2][1]) / s
	outAxis.Y = (m[2][0] - m[0][2]) / s
	outAxis.Z = (m[0][1] - m[1][0]) / s
}

// SetRotationEulers overwrites the current contents of this matrix with a rotation matrix
// that rotates first by yaw (y rotation), then by pitch (x rotation), and then by roll (z rotation)
//
// yawRad - Angle to rotate yaw in radians
//
// pitchRad - Angle to rotate pitch in radians
//
// rollRad - Angle to rotate roll in radians
func (m *Mat3x3[T]) SetRotationEulers(yawRad, pitchRad, rollRad float64) {
	yawCos := T(math.Cos(yawRad))
	pitchCos := T(math.Cos(pitchRad))
	rollCos := T(math.Cos(rollRad))
	yawSin := T(math.Sin(yawRad))
	pitchSin := T(math.Sin(pitchRad))
	rollSin := T(math.Sin(rollRad))

	m[0][0] = rollCos * yawCos
	m[0][1] = rollSin*yawCos + pitchSin*yawSin*rollCos
	m[0][2] = -yawSin * pitchCos
	m[1][0] = -rollSin * pitchCos
	m[1][1] = rollCos * pitchCos
	m[1][2] = pitchSin
	m[2][0] = yawSin
	m[2][1] = rollSin*yawSin - yawCos*pitchSin*rollCos
	m[2][2] = pitchCos * yawCos
}

// SetRotationY overwrites the current contents of this matrix with a rotation matrix that
// rotates around the y axis by the specified amount
//
// yawRad - The angle to rotate around the y axis in radians
func (m *Mat3x3[T]) SetRotationY(yawRad float64) {
	cos := T(math.Cos(yawRad))
	sin := T(math.Sin(yawRad))

	m00 := cos
	m02 := -sin

	m20 := sin
	m22 := cos

	m[0][0] = m00
	m[0][1] = 0
	m[0][2] = m02
	m[1][0] = 0
	m[1][1] = 1
	m[1][2] = 0
	m[2][0] = m20
	m[2][1] = 0
	m[2][2] = m22
}

// SetRotationX overwrites the current contents of this matrix with a rotation matrix that
// rotates around the x axis by the specified amount
//
// pitchRad - The angle to rotate around the x axis in radians
func (m *Mat3x3[T]) SetRotationX(pitchRad float64) {
	cos := T(math.Cos(pitchRad))
	sin := T(math.Sin(pitchRad))

	m11 := cos
	m12 := sin

	m21 := -sin
	m22 := cos

	m[0][0] = 1
	m[0][1] = 0
	m[0][2] = 0
	m[1][0] = 0
	m[1][1] = m11
	m[1][2] = m12
	m[2][0] = 0
	m[2][1] = m21
	m[2][2] = m22
}

// SetRotationZ overwrites the current contents of this matrix with a rotation matrix that
// rotates around the z axis by the specified amount
//
// rollRad - The angle to rotate around the z axis in radians
func (m *Mat3x3[T]) SetRotationZ(rollRad float64) {
	cos := T(math.Cos(rollRad))
	sin := T(math.Sin(rollRad))

	m00 := cos
	m01 := sin

	m10 := -sin
	m11 := cos

	m[0][0] = m00
	m[0][1] = m01
	m[0][2] = 0
	m[1][0] = m10
	m[1][1] = m11
	m[1][2] = 0
	m[2][0] = 0
	m[2][1] = 0
	m[2][2] = 1
}

// SetScale overwrites the current contents of this matrix with a transform matrix that scales
// vectors by the 3 provided scalars
//
// x - Factor to scale by along the x axis
//
// y - Factor to scale by along the y axis
//
// z - Factor to scale by along the z axis
func (m *Mat3x3[T]) SetScale(x, y, z T) {
	m[0][0] = x
	m[0][1] = 0
	m[0][2] = 0
	m[1][0] = 0
	m[1][1] = y
	m[1][2] = 0
	m[2][0] = 0
	m[2][1] = 0
	m[2][2] = z
}

// Transpose mirrors this matrix across the diagonal
func (m *Mat3x3[T]) Transpose() {
	m[1][0], m[0][1] = m[0][1], m[1][0]
	m[2][0], m[0][2] = m[0][2], m[2][0]
	m[2][1], m[1][2] = m[1][2], m[2][1]
}

// Inverse inverts this matrix. If this matrix has no valid inverse, some
// entries will be set to NaN
func (m *Mat3x3[T]) Inverse() {
	determinant := m[0][0]*(m[1][1]*m[2][2]-m[2][1]*m[1][2]) -
		m[1][0]*(m[0][1]*m[2][2]-m[2][1]*m[0][2]) +
		m[2][0]*(m[0][1]*m[1][2]-m[1][1]*m[0][2])
	oneOverDeterminant := 1.0 / determinant

	m00 := (m[1][1]*m[2][2] - m[2][1]*m[1][2]) * oneOverDeterminant
	m10 := -(m[1][0]*m[2][2] - m[2][0]*m[1][2]) * oneOverDeterminant
	m20 := (m[1][0]*m[2][1] - m[2][0]*m[1][1]) * oneOverDeterminant
	m01 := -(m[0][1]*m[2][2] - m[2][1]*m[0][2]) * oneOverDeterminant
	m11 := (m[0][0]*m[2][2] - m[2][0]*m[0][2]) * oneOverDeterminant
	m21 := -(m[0][0]*m[2][1] - m[2][0]*m[0][1]) * oneOverDeterminant
	m02 := (m[0][1]*m[1][2] - m[1][1]*m[0][2]) * oneOverDeterminant
	m12 := -(m[0][0]*m[1][2] - m[1][0]*m[0][2]) * oneOverDeterminant
	m22 := (m[0][0]*m[1][1] - m[1][0]*m[0][1]) * oneOverDeterminant

	m[0][0] = m00
	m[1][0] = m10
	m[2][0] = m20
	m[0][1] = m01
	m[1][1] = m11
	m[2][1] = m21
	m[0][2] = m02
	m[1][2] = m12
	m[2][2] = m22
}

// Orthonormalize this matrix. All 3 columns will be normalized and made to lie
// at right angles to one another. Only the first column is guaranteed to maintain
// its current direction.
func (m *Mat3x3[T]) Orthonormalize() {
	col0 := Vec3[T]{m[0][0], m[0][1], m[0][2]}
	col1 := Vec3[T]{m[1][0], m[1][1], m[1][2]}
	col2 := Vec3[T]{m[2][0], m[2][1], m[2][2]}

	col0.Normalize()
	dot0 := col0.DotProduct(&col1)

	col1.X -= col0.X * dot0
	col1.Y -= col0.Y * dot0
	col1.Z -= col0.Z * dot0

	col1.Normalize()

	dot1 := col1.DotProduct(&col2)
	dot0 = col0.DotProduct(&col2)

	col2.X -= col0.X*dot0 + col1.X*dot1
	col2.Y -= col0.Y*dot0 + col1.Y*dot1
	col2.Z -= col0.Z*dot0 + col1.Z*dot1
	col2.Normalize()

	m[0][0] = col0.X
	m[0][1] = col0.Y
	m[0][2] = col0.Z
	m[1][0] = col1.X
	m[1][1] = col1.Y
	m[1][2] = col1.Z
	m[2][0] = col2.X
	m[2][1] = col2.Y
	m[2][2] = col2.Z
}

// ApplyTransform applies the provided transform matrix to this transform matrix. "Apply
// Transform" is just a matrix multiply with the operands reversed.
//
// other - The transform matrix to apply
func (m *Mat3x3[T]) ApplyTransform(other *Mat3x3[T]) {
	m00 := other[0][0]*m[0][0] + other[1][0]*m[0][1] + other[2][0]*m[0][2]
	m10 := other[0][0]*m[1][0] + other[1][0]*m[1][1] + other[2][0]*m[1][2]
	m20 := other[0][0]*m[2][0] + other[1][0]*m[2][1] + other[2][0]*m[2][2]

	m01 := other[0][1]*m[0][0] + other[1][1]*m[0][1] + other[2][1]*m[0][2]
	m11 := other[0][1]*m[1][0] + other[1][1]*m[1][1] + other[2][1]*m[1][2]
	m21 := other[0][1]*m[2][0] + other[1][1]*m[2][1] + other[2][1]*m[2][2]

	m02 := other[0][2]*m[0][0] + other[1][2]*m[0][1] + other[2][2]*m[0][2]
	m12 := other[0][2]*m[1][0] + other[1][2]*m[1][1] + other[2][2]*m[1][2]
	m22 := other[0][2]*m[2][0] + other[1][2]*m[2][1] + other[2][2]*m[2][2]

	m[0][0] = m00
	m[0][1] = m01
	m[0][2] = m02
	m[1][0] = m10
	m[1][1] = m11
	m[1][2] = m12
	m[2][0] = m20
	m[2][1] = m21
	m[2][2] = m22
}

// MultMat3x3 multiplies this matrix with the provided matrix and updates this
// matrix with the results.  In Matrix Multiplication between transform matrices,
// the left matrix is unintuitively applied to the right matrix rather than the other way around.
// You may prefer to use ApplyTransform.
//
// other - The right hand side of the multiplication operation.
func (m *Mat3x3[T]) MultMat3x3(other *Mat3x3[T]) {
	m00 := m[0][0]*other[0][0] + m[1][0]*other[0][1] + m[2][0]*other[0][2]
	m10 := m[0][0]*other[1][0] + m[1][0]*other[1][1] + m[2][0]*other[1][2]
	m20 := m[0][0]*other[2][0] + m[1][0]*other[2][1] + m[2][0]*other[2][2]

	m01 := m[0][1]*other[0][0] + m[1][1]*other[0][1] + m[2][1]*other[0][2]
	m11 := m[0][1]*other[1][0] + m[1][1]*other[1][1] + m[2][1]*other[1][2]
	m21 := m[0][1]*other[2][0] + m[1][1]*other[2][1] + m[2][1]*other[2][2]

	m02 := m[0][2]*other[0][0] + m[1][2]*other[0][1] + m[2][2]*other[0][2]
	m12 := m[0][2]*other[1][0] + m[1][2]*other[1][1] + m[2][2]*other[1][2]
	m22 := m[0][2]*other[2][0] + m[1][2]*other[2][1] + m[2][2]*other[2][2]

	m[0][0] = m00
	m[0][1] = m01
	m[0][2] = m02
	m[1][0] = m10
	m[1][1] = m11
	m[1][2] = m12
	m[2][0] = m20
	m[2][1] = m21
	m[2][2] = m22
}

// Proj2D transforms the matrix as a planar projection matrix along the provided
// normal axis
//
// normal - The normal of the plane
func (m *Mat3x3[T]) Proj2D(normal *Vec3[T]) {
	transform00 := 1.0 - normal.X*normal.X
	transform01 := -normal.X * normal.Y
	transform10 := -normal.X * normal.Y
	transform11 := 1.0 - normal.Y*normal.Y

	m00 := m[0][0]*transform00 + m[0][1]*transform10
	m01 := m[0][0]*transform01 + m[0][1]*transform11

	m10 := m[1][0]*transform00 + m[1][1]*transform10
	m11 := m[1][0]*transform01 + m[1][1]*transform11

	m20 := m[2][0]*transform00 + m[2][1]*transform10
	m21 := m[2][0]*transform01 + m[2][1]*transform11

	m[0][0] = m00
	m[0][1] = m01
	m[1][0] = m10
	m[1][1] = m11
	m[2][0] = m20
	m[2][1] = m21
}

// SetShearX overwrites this matrix with a shear matrix that shears along the X
// axis by the provided factor scalars
//
// y - y shear factor
//
// z - z shear factor
func (m *Mat3x3[T]) SetShearX(y, z T) {
	m[0][0] = 1
	m[1][0] = 0
	m[2][0] = 0
	m[0][1] = y
	m[1][1] = 1
	m[2][1] = 0
	m[0][2] = z
	m[1][2] = 0
	m[2][2] = 1
}

// ShearX applies a shear transform to this matrix along the X axis by the provided
// factor scalars
//
// y - y shear factor
//
// z - z shear factor
func (m *Mat3x3[T]) ShearX(y, z T) {

	m01 := y*m[0][0] + m[0][1]
	m11 := y*m[1][0] + m[1][1]
	m21 := y*m[2][0] + m[2][1]

	m02 := z*m[0][0] + m[0][2]
	m12 := z*m[1][0] + m[1][2]
	m22 := z*m[2][0] + m[2][2]

	m[0][1] = m01
	m[0][2] = m02
	m[1][1] = m11
	m[1][2] = m12
	m[2][1] = m21
	m[2][2] = m22
}

// SetShearY overwrites this matrix with a shear matrix that shears along the Y
// axis by the provided factor scalars
//
// x - x shear factor
//
// z - z shear factor
func (m *Mat3x3[T]) SetShearY(x, z T) {
	m[0][0] = 1
	m[1][0] = x
	m[2][0] = 0
	m[0][1] = 0
	m[1][1] = 1
	m[2][1] = 0
	m[0][2] = 0
	m[1][2] = z
	m[2][2] = 1
}

// ShearY applies a shear transform to this matrix along the Y axis by the provided
// factor scalars
//
// x - x shear factor
//
// z - z shear factor
func (m *Mat3x3[T]) ShearY(x, z T) {
	m00 := m[0][0] + x*m[0][1]
	m10 := m[1][0] + x*m[1][1]
	m20 := m[2][0] + x*m[2][1]

	m02 := z*m[0][1] + m[0][2]
	m12 := z*m[1][1] + m[1][2]
	m22 := z*m[2][1] + m[2][2]

	m[0][0] = m00
	m[0][2] = m02
	m[1][0] = m10
	m[1][2] = m12
	m[2][0] = m20
	m[2][2] = m22
}

// SetShearZ overwrites this matrix with a shear matrix that shears along the Z
// axis by the provided factor scalars
//
// x - x shear factor
//
// y - y shear factor
func (m *Mat3x3[T]) SetShearZ(x, y T) {
	m[0][0] = 1
	m[1][0] = 0
	m[2][0] = x
	m[0][1] = 0
	m[1][1] = 1
	m[2][1] = y
	m[0][2] = 0
	m[1][2] = 0
	m[2][2] = 1
}

// ShearZ applies a shear transform to this matrix along the Z axis by the provided
// factor scalars
//
// x - x shear factor
//
// y - y shear factor
func (m *Mat3x3[T]) ShearZ(x, y T) {

	m00 := m[0][0] + x*m[0][2]
	m10 := m[1][0] + x*m[1][2]
	m20 := m[2][0] + x*m[2][2]

	m01 := m[0][1] + y*m[0][2]
	m11 := m[1][1] + y*m[1][2]
	m21 := m[2][1] + y*m[2][2]

	m[0][0] = m00
	m[0][1] = m01
	m[1][0] = m10
	m[1][1] = m11
	m[2][0] = m20
	m[2][1] = m21
}

// Scale applies a scale transform to this matrix that scales vectors by the 3 provided scalars
//
// x - Factor to scale by along the x axis
//
// y - Factor to scale by along the y axis
//
// z - Factor to scale by along the z axis
func (m *Mat3x3[T]) Scale(x, y, z T) {
	m[0][0] *= x
	m[0][1] *= y
	m[0][2] *= z
	m[1][0] *= x
	m[1][1] *= y
	m[1][2] *= z
	m[2][0] *= x
	m[2][1] *= y
	m[2][2] *= z
}

// SetRotationAroundAxis overwrites the current contents of this matrix with a rotation matrix that rotates
// around the provided axis by the provided angle in radians.
//
// axis - A 3-element vector that is normal to the angle of rotation. It does not need to be normalized.
//
// angleRad - The amount to rotate in radians
func (m *Mat3x3[T]) SetRotationAroundAxis(axis *Vec3[T], angleRad float64) {
	var unitAxis Vec3[T]
	unitAxis.SetVec3(axis)
	unitAxis.Normalize()

	cos := T(math.Cos(angleRad))
	sin := T(math.Sin(angleRad))

	inverseCos := 1 - cos

	m[0][0] = cos + unitAxis.X*unitAxis.X*inverseCos
	m[0][1] = unitAxis.X*unitAxis.Y*inverseCos + sin*unitAxis.Z
	m[0][2] = unitAxis.X*unitAxis.Z*inverseCos - sin*unitAxis.Y
	m[1][0] = unitAxis.Y*unitAxis.X*inverseCos - sin*unitAxis.Z
	m[1][1] = cos + unitAxis.Y*unitAxis.Y*inverseCos
	m[1][2] = unitAxis.Y*unitAxis.Z*inverseCos + sin*unitAxis.X
	m[2][0] = unitAxis.Z*unitAxis.X*inverseCos + sin*unitAxis.Y
	m[2][1] = unitAxis.Z*unitAxis.Y*inverseCos - sin*unitAxis.X
	m[2][2] = cos + unitAxis.Z*unitAxis.Z*inverseCos
}

// RotateAroundAxis transforms this matrix by applying a rotation around the provided axis by the
// provided angle in radians
//
// axis - A 3-element vector that is normal to the angle of rotation. It does not need to be normalized.
//
// angleRad - The amount to rotate in radians
func (m *Mat3x3[T]) RotateAroundAxis(axis *Vec3[T], angleRad float64) {
	var unitAxis Vec3[T]
	unitAxis.SetVec3(axis)
	unitAxis.Normalize()

	cos := T(math.Cos(angleRad))
	sin := T(math.Sin(angleRad))

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

	m00 := rotate00*m[0][0] + rotate10*m[0][1] + rotate20*m[0][2]
	m01 := rotate01*m[0][0] + rotate11*m[0][1] + rotate21*m[0][2]
	m02 := rotate02*m[0][0] + rotate12*m[0][1] + rotate22*m[0][2]

	m10 := rotate00*m[1][0] + rotate10*m[1][1] + rotate20*m[1][2]
	m11 := rotate01*m[1][0] + rotate11*m[1][1] + rotate21*m[1][2]
	m12 := rotate02*m[1][0] + rotate12*m[1][1] + rotate22*m[1][2]

	m20 := rotate00*m[2][0] + rotate10*m[2][1] + rotate20*m[2][2]
	m21 := rotate01*m[2][0] + rotate11*m[2][1] + rotate21*m[2][2]
	m22 := rotate02*m[2][0] + rotate12*m[2][1] + rotate22*m[2][2]

	m[0][0] = m00
	m[0][1] = m01
	m[0][2] = m02
	m[1][0] = m10
	m[1][1] = m11
	m[1][2] = m12
	m[2][0] = m20
	m[2][1] = m21
	m[2][2] = m22
}

// RotateX applies a transformation to this matrix that rotates around the x axis by the specified amount
//
// angleRad - The angle to rotate around the x axis in radians
func (m *Mat3x3[T]) RotateX(angleRad float64) {
	cos := T(math.Cos(angleRad))
	sin := T(math.Sin(angleRad))

	m01 := cos*m[0][1] - sin*m[0][2]
	m11 := cos*m[1][1] - sin*m[1][2]
	m21 := cos*m[2][1] - sin*m[2][2]

	m02 := sin*m[0][1] + cos*m[0][2]
	m12 := sin*m[1][1] + cos*m[1][2]
	m22 := sin*m[2][1] + cos*m[2][2]

	m[0][1] = m01
	m[1][1] = m11
	m[2][1] = m21
	m[0][2] = m02
	m[1][2] = m12
	m[2][2] = m22
}

// RotateY applies a transformation to this matrix that rotates around the y axis by the specified amount
//
// angleRad - The angle to rotate around the y axis in radians
func (m *Mat3x3[T]) RotateY(angleRad float64) {
	cos := T(math.Cos(angleRad))
	sin := T(math.Sin(angleRad))

	m00 := cos*m[0][0] + sin*m[0][2]
	m10 := cos*m[1][0] + sin*m[1][2]
	m20 := cos*m[2][0] + sin*m[2][2]

	m02 := -sin*m[0][0] + cos*m[0][2]
	m12 := -sin*m[1][0] + cos*m[1][2]
	m22 := -sin*m[2][0] + cos*m[2][2]

	m[0][0] = m00
	m[1][0] = m10
	m[2][0] = m20
	m[0][2] = m02
	m[1][2] = m12
	m[2][2] = m22
}

// RotateZ applies a transformation to this matrix that rotates around the z axis by the specified amount
//
// angleRad - The angle to rotate around the z axis in radians
func (m *Mat3x3[T]) RotateZ(angleRad float64) {
	cos := T(math.Cos(angleRad))
	sin := T(math.Sin(angleRad))

	m00 := cos*m[0][0] - sin*m[0][1]
	m10 := cos*m[1][0] - sin*m[1][1]
	m20 := cos*m[2][0] - sin*m[2][1]

	m01 := sin*m[0][0] + cos*m[0][1]
	m11 := sin*m[1][0] + cos*m[1][1]
	m21 := sin*m[2][0] + cos*m[2][1]

	m[0][0] = m00
	m[1][0] = m10
	m[2][0] = m20
	m[0][1] = m01
	m[1][1] = m11
	m[2][1] = m21
}

// IsNormalized returns true if every row and every column in this matrix has a length of 1
//
// epsilon - The epsilon value to use in floating point comparisons. This much floating point
// drift is permitted before the method returns false. 0.0001 is a common epsilon value
func (m *Mat3x3[T]) IsNormalized(epsilon T) bool {
	for i := 0; i < 3; i++ {
		column := abs[T](m[i][0]*m[i][0] + m[i][1]*m[i][1] + m[i][2]*m[i][2])
		if column-T(1) > T(2)*epsilon {
			return false
		}
	}

	for i := 0; i < 3; i++ {
		row := abs[T](m[0][i]*m[0][i] + m[1][i]*m[1][i] + m[2][i]*m[2][i])
		for row-T(1) > T(2)*epsilon {
			return false
		}
	}

	return true
}

// IsNull returns true if every column in this matrix has a length of 1
//
// epsilon - The epsilon value to use in floating point comparisons. This much floating point
// drift is permitted before the method returns false. 0.0001 is a common epsilon value
func (m *Mat3x3[T]) IsNull(epsilon T) bool {
	for i := 0; i < 3; i++ {
		column := abs[T](m[i][0]*m[i][0] + m[i][1]*m[i][1] + m[i][2]*m[i][2])
		if column > epsilon {
			return false
		}
	}

	return true
}

// Equal returns true if every entry in this matrix is equal to every entry in the provided matrix
//
// other - The matrix to compare this matrix to
//
// epsilon - The epsilon value to use in floating point comparisons. This much floating point
// drift is permitted before the method returns false. 0.0001 is a common epsilon value
func (m *Mat3x3[T]) Equal(other *Mat3x3[T], epsilon T) bool {
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			if abs[T](m[i][j]-other[i][j]) > epsilon {
				return false
			}
		}
	}

	return true
}
