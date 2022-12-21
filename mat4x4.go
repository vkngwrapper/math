package math

import "math"

type Mat4Row[T FloatingPoint] [4]T

// Mat4x4 is a two-dimensional array of floating point-compatible values that can be used
// for 4x4 matrix arithmetic and all 3d matrix transformations
type Mat4x4[T FloatingPoint] [4]Mat4Row[T]

// SetIdentity overwrites the current contents of the matrix with the identity matrix
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

// SetDiagonalScalar overwrites the current contents of the matrix with an identity matrix,
// multiplied by the provided scalar value
//
// s - The scalar value to multiply the matrix entries by
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

// SetDiagonalVector overwrites the current contents of the matrix with an identity matrix,
// in which the 1 values are replaced with the contents of a provided vector, with each vector
// element corresponding to a column in the matrix
//
// v - The vector containing the diagonal entries that will be populated into the matrix
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

// SetFrustum overwrites the current matrix values with a projection matrix for a frustum
// with the provided properties
//
// left - The left boundary of the view frustum
// right - The right boundary of the view frustum
// bottom - The bottom boundary of the view frustum
// top - The top boundary of the view frustum
// nearVal - The near clipping plane distance from origin
// farVal - The far clipping plane distance from origin
func (m *Mat4x4[T]) SetFrustum(left, right, bottom, top, nearVal, farVal T) *Mat4x4[T] {
	m[0][0] = (2 * nearVal) / (right - left)
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

// SetOrientation overwrites the current matrix values with a rotation matrix which rotates
// from the provided source vector to the provided target vector
//
// source - A unit vector indicating the starting direction of the rotation
//
// target - A unit vector indicating the target direction of the rotation
func (m *Mat4x4[T]) SetOrientation(source *Vec3[T], target *Vec3[T]) *Mat4x4[T] {
	if target.Equal(source, 0.0001) {
		m.SetIdentity()
		return m
	}

	var rotationAxis Vec3[T]
	rotationAxis.SetVec3(source).CrossProduct(target)

	angle := math.Acos(float64(target.DotProduct(source)))

	return m.SetRotateAroundAxis(&rotationAxis, angle)
}

// SetTranslation overwrites the current matrix values with a transform matrix which translates
// vectors by the provided values
//
// x - The amount to translate along the x axis
//
// y - The amount to translate along the y axis
//
// z - The amount to translate along the z axis
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

// SetScale overwrites the current contents of this matrix with a transform matrix that scales
// vectors by the 3 provided scalars
//
// x - Factor to scale by along the x axis
//
// y - Factor to scale by along the y axis
//
// z - Factor to scale by along the z axis
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

// SetMat4x4 overwrites the current contents of the matrix with the contents of the provided
// matrix
//
// other - The matrix to initialize from
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

// SetMat3x3 overwrites the current contents of the matrix with the identity matrix and
// layers the provided 3x3 matrix over the upper left 3x3 area
//
// other - The 3x3 matrix to initialize from
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

// SetMat2x2 overwrites the current contents of the matrix with the identity matrix and
// layers the provided 2x2 matrix over the upper left 2x2 area of the matrix
//
// other - The 2x2 matrix to initialize from
func (m *Mat4x4[T]) SetMat2x2(other *Mat2x2[T]) *Mat4x4[T] {
	m[0][0] = other[0][0]
	m[0][1] = other[0][1]
	m[0][2] = 0
	m[0][3] = 0
	m[1][0] = other[1][0]
	m[1][1] = other[1][1]
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

// SetQuaternion overwrites the current contents of the matrix with an affine rotation matrix
// that contains the rotation of the provided Quaternion
//
// other - The quaternion to initialize from
func (m *Mat4x4[T]) SetQuaternion(other *Quaternion[T]) *Mat4x4[T] {
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
	m[0][3] = 0
	m[1][0] = 2.0 * (xy - wz)
	m[1][1] = 1.0 - 2.0*(xx+zz)
	m[1][2] = 2.0 * (yz + wx)
	m[1][3] = 0
	m[2][0] = 2.0 * (xz + wy)
	m[2][1] = 2.0 * (yz - wx)
	m[2][2] = 1.0 - 2.0*(xx+yy)
	m[2][3] = 0
	m[3][0] = 0
	m[3][1] = 0
	m[3][2] = 0
	m[3][3] = 1

	return m
}

// SetColMajor overwrites the current contents of the matrix with 4 matrix columns
// passed in as 4-element vectors
//
// c1 - The first column of the matrix
// c2 - The second column of the matrix
// c3 - The third column of the matrix
// c4 - The fourth column of the matrix
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

// SetRowMajor overwrites the current contents of the matrix with 4 matrix rows
// passed in as 4-element vectors
//
// r1 - The first row of the matrix
// r2 - The second row of the matrix
// r3 - The third row of the matrix
// r4 - The fourth row of the matrix
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

// SetRotationEulers overwrites the current contents of this matrix with a rotation matrix
// that rotates first by yaw (y rotation), then by pitch (x rotation), and then by roll (z rotation)
//
// yawRad - Angle to rotate yaw in radians
//
// pitchRad - Angle to rotate pitch in radians
//
// rollRad - Angle to rotate roll in radians
func (m *Mat4x4[T]) SetRotationEulers(yawRad, pitchRad, rollRad float64) *Mat4x4[T] {
	yawCos := T(math.Cos(yawRad))
	pitchCos := T(math.Cos(pitchRad))
	rollCos := T(math.Cos(rollRad))
	yawSin := T(math.Sin(yawRad))
	pitchSin := T(math.Sin(pitchRad))
	rollSin := T(math.Sin(rollRad))

	m[0][0] = rollCos * yawCos
	m[0][1] = rollSin*yawCos + pitchSin*yawSin*rollCos
	m[0][2] = -yawSin * pitchCos
	m[0][3] = 0
	m[1][0] = -rollSin * pitchCos
	m[1][1] = rollCos * pitchCos
	m[1][2] = pitchSin
	m[1][3] = 0
	m[2][0] = yawSin
	m[2][1] = rollSin*yawSin - yawCos*pitchSin*rollCos
	m[2][2] = pitchCos * yawCos
	m[2][3] = 0
	m[3][0] = 0
	m[3][1] = 0
	m[3][2] = 0
	m[3][3] = 1

	return m
}

// SetRotationY overwrites the current contents of this matrix with a rotation matrix that
// rotates around the y axis by the specified amount
//
// yawRad - The angle to rotate around the y axis in radians
func (m *Mat4x4[T]) SetRotationY(yawRad float64) *Mat4x4[T] {
	cos := T(math.Cos(yawRad))
	sin := T(math.Sin(yawRad))

	m00 := cos
	m02 := -sin

	m20 := sin
	m22 := cos

	m[0][0] = m00
	m[0][1] = 0
	m[0][2] = m02
	m[0][3] = 0
	m[1][0] = 0
	m[1][1] = 1
	m[1][2] = 0
	m[1][3] = 0
	m[2][0] = m20
	m[2][1] = 0
	m[2][2] = m22
	m[2][3] = 0
	m[3][0] = 0
	m[3][1] = 0
	m[3][2] = 0
	m[3][3] = 1

	return m
}

// SetRotationX overwrites the current contents of this matrix with a rotation matrix that
// rotates around the x axis by the specified amount
//
// pitchRad - The angle to rotate around the x axis in radians
func (m *Mat4x4[T]) SetRotationX(pitchRad float64) *Mat4x4[T] {
	cos := T(math.Cos(pitchRad))
	sin := T(math.Sin(pitchRad))

	m11 := cos
	m12 := sin

	m21 := -sin
	m22 := cos

	m[0][0] = 1
	m[0][1] = 0
	m[0][2] = 0
	m[0][3] = 0
	m[1][0] = 0
	m[1][1] = m11
	m[1][2] = m12
	m[1][3] = 0
	m[2][0] = 0
	m[2][1] = m21
	m[2][2] = m22
	m[2][3] = 0
	m[3][0] = 0
	m[3][1] = 0
	m[3][2] = 0
	m[3][3] = 1

	return m
}

// SetRotationZ overwrites the current contents of this matrix with a rotation matrix that
// rotates around the z axis by the specified amount
//
// rollRad - The angle to rotate around the z axis in radians
func (m *Mat4x4[T]) SetRotationZ(rollRad float64) *Mat4x4[T] {
	cos := T(math.Cos(rollRad))
	sin := T(math.Sin(rollRad))

	m00 := cos
	m01 := sin

	m10 := -sin
	m11 := cos

	m[0][0] = m00
	m[0][1] = m01
	m[0][2] = 0
	m[0][3] = 0
	m[1][0] = m10
	m[1][1] = m11
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

// SetAffineMat4x4 overwrites the current contents of this matrix with the identity matrix,
// and then copies the upper left 3x3 portion of the provided 4x4 matrix into the upper left
// 3x3 portion of this matrix.
//
// other - The 4x4 matrix to copy from
func (m *Mat4x4[T]) SetAffineMat4x4(other *Mat4x4[T]) *Mat4x4[T] {
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

// SetInfinitePerspective overwrites the current contents of this matrix with a projection matrix
// that is used to render objects that are infinitely-far away
//
// fovYRad - The vertical field of view of the projection in radians
//
// aspectRatio - The aspect ratio of the viewport
//
// zNear - The near clipping plane's distance from the origin
func (m *Mat4x4[T]) SetInfinitePerspective(fovYRad float64, aspectRatio, zNear T) *Mat4x4[T] {
	r := T(math.Tan(fovYRad/2)) * zNear
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
	m[2][2] = -0.9999
	m[2][3] = -1
	m[3][0] = 0
	m[3][1] = 0
	m[3][2] = -1.9999 * zNear
	m[3][3] = 0

	return m
}

// SetOrthographic2D overwrites the current contents of this matrix with a projection matrix that is
// used to render 2-dimensional objects with no perspective effect and no depth
//
// left - The left boundary of the viewport
//
// right - The right boundary of the viewport
//
// bottom - The bottom boundary of the viewport
//
// top - The top boundary of the viewport
func (m *Mat4x4[T]) SetOrthographic2D(left, right, bottom, top T) *Mat4x4[T] {
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

// SetOrthographic overwrites the current contents of this matrix with a projection matrix that is used
// to render objects with no perspective, generally used for 2d rendering
//
// left - The left boundary of the viewport
//
// right - The right boundary of the viewport
//
// bottom - The bottom boundary of the viewport
//
// top - The top boundary of the viewport
//
// zNear - The near clipping plane's distance from the origin
//
// zFar - The far clipping plane's distance from the origin
func (m *Mat4x4[T]) SetOrthographic(left, right, bottom, top, zNear, zFar T) *Mat4x4[T] {
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

// SetPerspective overwrites the current contents of this matrix with a projection matrix that
// is used to render in a 3d perspective
//
// fovYRad - The vertical field of view of the projection, in radians
//
// aspectRatio - The aspect ratio of the viewport
//
// zNear - The near clipping plane's distance from the origin
//
// zFar - The far clipping plane's distance from the origin
func (m *Mat4x4[T]) SetPerspective(fovYRad float64, aspectRatio, zNear, zFar T) *Mat4x4[T] {
	tanHalfFovY := T(math.Tan(fovYRad * 0.5))

	m[0][0] = 1 / (aspectRatio * tanHalfFovY)
	m[0][1] = 0
	m[0][2] = 0
	m[0][3] = 0
	m[1][0] = 0
	m[1][1] = -1 / tanHalfFovY
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

// SetLookAt overwrites the current contents of this matrix with a view matrix that is used
// to view a scene, given a particular camera position and orientation
//
// eyePosition - A 3d vector indicating the location of the camera
//
// target - A 3d vector indicating the point the camera is focusing on
//
// up - A 3d vector indicating the camera's up direction
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

// SetMultMat4x4 multiplies two 4x4 matrices together and overwrites the current contents
// of this matrix with the results. In Matrix Multiplication between transform matrices,
// the left matrix is unintuitively applied to the right matrix rather than the other way around.
// You may prefer to use ApplyTransform.
//
// lhs - The left operand of the multiplication operation
// rhs - The right operand of the multiplication operation
func (m *Mat4x4[T]) SetMultMat4x4(lhs, rhs *Mat4x4[T]) *Mat4x4[T] {
	m[0][0] = lhs[0][0]*rhs[0][0] + lhs[1][0]*rhs[0][1] + lhs[2][0]*rhs[0][2] + lhs[3][0]*rhs[0][3]
	m[0][1] = lhs[0][1]*rhs[0][0] + lhs[1][1]*rhs[0][1] + lhs[2][1]*rhs[0][2] + lhs[3][1]*rhs[0][3]
	m[0][2] = lhs[0][2]*rhs[0][0] + lhs[1][2]*rhs[0][1] + lhs[2][2]*rhs[0][2] + lhs[3][2]*rhs[0][3]
	m[0][3] = lhs[0][3]*rhs[0][0] + lhs[1][3]*rhs[0][1] + lhs[2][3]*rhs[0][2] + lhs[3][3]*rhs[0][3]
	m[1][0] = lhs[0][0]*rhs[1][0] + lhs[1][0]*rhs[1][1] + lhs[2][0]*rhs[1][2] + lhs[3][0]*rhs[1][3]
	m[1][1] = lhs[0][1]*rhs[1][0] + lhs[1][1]*rhs[1][1] + lhs[2][1]*rhs[1][2] + lhs[3][1]*rhs[1][3]
	m[1][2] = lhs[0][2]*rhs[1][0] + lhs[1][2]*rhs[1][1] + lhs[2][2]*rhs[1][2] + lhs[3][2]*rhs[1][3]
	m[1][3] = lhs[0][3]*rhs[1][0] + lhs[1][3]*rhs[1][1] + lhs[2][3]*rhs[1][2] + lhs[3][3]*rhs[1][3]
	m[2][0] = lhs[0][0]*rhs[2][0] + lhs[1][0]*rhs[2][1] + lhs[2][0]*rhs[2][2] + lhs[3][0]*rhs[2][3]
	m[2][1] = lhs[0][1]*rhs[2][0] + lhs[1][1]*rhs[2][1] + lhs[2][1]*rhs[2][2] + lhs[3][1]*rhs[2][3]
	m[2][2] = lhs[0][2]*rhs[2][0] + lhs[1][2]*rhs[2][1] + lhs[2][2]*rhs[2][2] + lhs[3][2]*rhs[2][3]
	m[2][3] = lhs[0][3]*rhs[2][0] + lhs[1][3]*rhs[2][1] + lhs[2][3]*rhs[2][2] + lhs[3][3]*rhs[2][3]
	m[3][0] = lhs[0][0]*rhs[3][0] + lhs[1][0]*rhs[3][1] + lhs[2][0]*rhs[3][2] + lhs[3][0]*rhs[3][3]
	m[3][1] = lhs[0][1]*rhs[3][0] + lhs[1][1]*rhs[3][1] + lhs[2][1]*rhs[3][2] + lhs[3][1]*rhs[3][3]
	m[3][2] = lhs[0][2]*rhs[3][0] + lhs[1][2]*rhs[3][1] + lhs[2][2]*rhs[3][2] + lhs[3][2]*rhs[3][3]
	m[3][3] = lhs[0][3]*rhs[3][0] + lhs[1][3]*rhs[3][1] + lhs[2][3]*rhs[3][2] + lhs[3][3]*rhs[3][3]

	return m
}

// SetInterpolateMat4x4 overwrites the current contents of this matrix with a transform
// matrix created by interpolating between the rotation and translation of two transform
// matrices.
//
// lhs - The "start" transform matrix in this interpolation
//
// rhs - The "end" transform matrix in this interpolation
//
// delta - A value between 0 and 1 indicating how far to interpolate between the two matrices
func (m *Mat4x4[T]) SetInterpolateMat4x4(lhs, rhs *Mat4x4[T], delta T) *Mat4x4[T] {
	var oldMatrix Mat4x4[T]
	oldMatrix.SetMat4x4(lhs)

	var thisRotation Mat4x4[T]
	thisRotation.SetAffineMat4x4(lhs)

	var transposedRotation Mat4x4[T]
	transposedRotation.SetMat4x4(&thisRotation).Transpose()

	var deltaRotation Mat4x4[T]
	deltaRotation.SetMultMat4x4(rhs, &transposedRotation)

	var deltaAxis Vec3[T]
	var deltaAngle float64
	deltaRotation.GetAxisAngle(&deltaAxis, &deltaAngle)

	m.SetRotateAroundAxis(&deltaAxis, deltaAngle).MultMat4x4(&thisRotation)
	m[3][0] = oldMatrix[3][0] + delta*(rhs[3][0]-oldMatrix[3][0])
	m[3][1] = oldMatrix[3][1] + delta*(rhs[3][1]-oldMatrix[3][1])
	m[3][2] = oldMatrix[3][2] + delta*(rhs[3][2]-oldMatrix[3][2])

	return m
}

// GetAxisAngle retrieves the axis and angle of rotation for this matrix, if it is a rotation matrix.
// If it is not a rotation matrix, the results are undefined.
//
// outAxis - A pointer to a 3-element vector that will be populated with a unit vector normal
// to the angle of rotation.
//
// outAngleRad - A pointer to a float64 that will be populated with the amount to rotate in radians
func (m *Mat4x4[T]) GetAxisAngle(outAxis *Vec3[T], outAngleRad *float64) {
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

// Transpose mirrors this matrix across the diagonal
func (m *Mat4x4[T]) Transpose() *Mat4x4[T] {
	m[1][0], m[0][1] = m[0][1], m[1][0]
	m[2][0], m[0][2] = m[0][2], m[2][0]
	m[2][1], m[1][2] = m[1][2], m[2][1]

	m[3][0], m[0][3] = m[0][3], m[3][0]
	m[3][1], m[1][3] = m[1][3], m[3][1]
	m[3][2], m[2][3] = m[2][3], m[3][2]

	return m
}

// Inverse inverts this matrix. If this matrix has no valid inverse, some
// entries will be set to NaN
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

	m00 := m[1][1]*coefficient0 - m[1][2]*coefficient4 + m[1][3]*coefficient8
	m01 := -m[0][1]*coefficient0 + m[0][2]*coefficient4 - m[0][3]*coefficient8
	m02 := m[0][1]*coefficient2 - m[0][2]*coefficient6 + m[0][3]*coefficient10
	m03 := -m[0][1]*coefficient3 + m[0][2]*coefficient7 - m[0][3]*coefficient11
	m10 := -m[1][0]*coefficient0 + m[1][2]*coefficient12 - m[1][3]*coefficient16
	m11 := m[0][0]*coefficient0 - m[0][2]*coefficient12 + m[0][3]*coefficient16
	m12 := -m[0][0]*coefficient2 + m[0][2]*coefficient14 - m[0][3]*coefficient18
	m13 := m[0][0]*coefficient3 - m[0][2]*coefficient15 + m[0][3]*coefficient19
	m20 := m[1][0]*coefficient4 - m[1][1]*coefficient12 + m[1][3]*coefficient20
	m21 := -m[0][0]*coefficient4 + m[0][1]*coefficient12 - m[0][3]*coefficient20
	m22 := m[0][0]*coefficient6 - m[0][1]*coefficient14 + m[0][3]*coefficient22
	m23 := -m[0][0]*coefficient7 + m[0][1]*coefficient15 - m[0][3]*coefficient23
	m30 := -m[1][0]*coefficient8 + m[1][1]*coefficient16 - m[1][2]*coefficient20
	m31 := m[0][0]*coefficient8 - m[0][1]*coefficient16 + m[0][2]*coefficient20
	m32 := -m[0][0]*coefficient10 + m[0][1]*coefficient18 - m[0][2]*coefficient22
	m33 := m[0][0]*coefficient11 - m[0][1]*coefficient19 + m[0][2]*coefficient23

	determinant := m[0][0]*m00 + m[0][1]*m10 + m[0][2]*m20 + m[0][3]*m30

	var oneOverDeterminant T
	if determinant < 0.0001 {
		oneOverDeterminant = T(math.NaN())
	} else {
		oneOverDeterminant = 1.0 / determinant
	}

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

// MultMat4x4 multiplies this matrix with the provided matrix and updates this
// matrix with the results.  In Matrix Multiplication between transform matrices,
// the left matrix is unintuitively applied to the right matrix rather than the other way around.
// You may prefer to use ApplyTransform.
//
// other - The right hand side of the multiplication operation.
func (m *Mat4x4[T]) MultMat4x4(other *Mat4x4[T]) *Mat4x4[T] {
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

// ApplyTransform applies the provided transform matrix to this transform matrix. "Apply
// Transform" is just a matrix multiply with the operands reversed.
//
// other - The transform matrix to apply
func (m *Mat4x4[T]) ApplyTransform(transform *Mat4x4[T]) *Mat4x4[T] {
	m00 := transform[0][0]*m[0][0] + transform[1][0]*m[0][1] + transform[2][0]*m[0][2] + transform[3][0]*m[0][3]
	m10 := transform[0][0]*m[1][0] + transform[1][0]*m[1][1] + transform[2][0]*m[1][2] + transform[3][0]*m[1][3]
	m20 := transform[0][0]*m[2][0] + transform[1][0]*m[2][1] + transform[2][0]*m[2][2] + transform[3][0]*m[2][3]
	m30 := transform[0][0]*m[3][0] + transform[1][0]*m[3][1] + transform[2][0]*m[3][2] + transform[3][0]*m[3][3]

	m01 := transform[0][1]*m[0][0] + transform[1][1]*m[0][1] + transform[2][1]*m[0][2] + transform[3][1]*m[0][3]
	m11 := transform[0][1]*m[1][0] + transform[1][1]*m[1][1] + transform[2][1]*m[1][2] + transform[3][1]*m[1][3]
	m21 := transform[0][1]*m[2][0] + transform[1][1]*m[2][1] + transform[2][1]*m[2][2] + transform[3][1]*m[2][3]
	m31 := transform[0][1]*m[3][0] + transform[1][1]*m[3][1] + transform[2][1]*m[3][2] + transform[3][1]*m[3][3]

	m02 := transform[0][2]*m[0][0] + transform[1][2]*m[0][1] + transform[2][2]*m[0][2] + transform[3][2]*m[0][3]
	m12 := transform[0][2]*m[1][0] + transform[1][2]*m[1][1] + transform[2][2]*m[1][2] + transform[3][2]*m[1][3]
	m22 := transform[0][2]*m[2][0] + transform[1][2]*m[2][1] + transform[2][2]*m[2][2] + transform[3][2]*m[2][3]
	m32 := transform[0][2]*m[3][0] + transform[1][2]*m[3][1] + transform[2][2]*m[3][2] + transform[3][2]*m[3][3]

	m03 := transform[0][3]*m[0][0] + transform[1][3]*m[0][1] + transform[2][3]*m[0][2] + transform[3][3]*m[0][3]
	m13 := transform[0][3]*m[1][0] + transform[1][3]*m[1][1] + transform[2][3]*m[1][2] + transform[3][3]*m[1][3]
	m23 := transform[0][3]*m[2][0] + transform[1][3]*m[2][1] + transform[2][3]*m[2][2] + transform[3][3]*m[2][3]
	m33 := transform[0][3]*m[3][0] + transform[1][3]*m[3][1] + transform[2][3]*m[3][2] + transform[3][3]*m[3][3]

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

// Proj3D transforms the matrix as a planar projection matrix along the provided
// normal axis
//
// normal - The normal of the plane
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

	return m.MultMat4x4(&proj)
}

// Translate applies a translation transform to this matrix that translates vectors
// by the provided scalar values
//
// x - Amount to translate along the x axis
//
// y - Amount to translate along the y axis
//
// z - Amount to translate along the z axis
func (m *Mat4x4[T]) Translate(x, y, z T) *Mat4x4[T] {
	m00 := m[0][0] + x*m[0][3]
	m10 := m[1][0] + x*m[1][3]
	m20 := m[2][0] + x*m[2][3]
	m30 := m[3][0] + x*m[3][3]

	m01 := m[0][1] + y*m[0][3]
	m11 := m[1][1] + y*m[1][3]
	m21 := m[2][1] + y*m[2][3]
	m31 := m[3][1] + y*m[3][3]

	m02 := m[0][2] + z*m[0][3]
	m12 := m[1][2] + z*m[1][3]
	m22 := m[2][2] + z*m[2][3]
	m32 := m[3][2] + z*m[3][3]

	m[0][0] = m00
	m[1][0] = m10
	m[2][0] = m20
	m[3][0] = m30
	m[0][1] = m01
	m[1][1] = m11
	m[2][1] = m21
	m[3][1] = m31
	m[0][2] = m02
	m[1][2] = m12
	m[2][2] = m22
	m[3][2] = m32

	return m
}

// Scale applies a scale transform to this matrix that scales vectors by the 3 provided scalars
//
// x - Factor to scale by along the x axis
//
// y - Factor to scale by along the y axis
//
// z - Factor to scale by along the z axis
func (m *Mat4x4[T]) Scale(x, y, z T) *Mat4x4[T] {
	m[0][0] *= x
	m[0][1] *= y
	m[0][2] *= z
	m[1][0] *= x
	m[1][1] *= y
	m[1][2] *= z
	m[2][0] *= x
	m[2][1] *= y
	m[2][2] *= z
	m[3][0] *= x
	m[3][1] *= y
	m[3][2] *= z

	return m
}

// SetRotateAroundAxis overwrites the current contents of this matrix with a rotation matrix that rotates
// around the provided axis by the provided angle in radians.
//
// axis - A 3-element vector that is normal to the angle of rotation. It does not need to be normalized.
//
// angleRad - The amount to rotate in radians
func (m *Mat4x4[T]) SetRotateAroundAxis(axis *Vec3[T], angleRad float64) *Mat4x4[T] {
	var unitAxis Vec3[T]
	unitAxis.SetVec3(axis).Normalize()

	cos := T(math.Cos(angleRad))
	sin := T(math.Sin(angleRad))

	inverseCos := 1 - cos

	m[0][0] = cos + unitAxis.X*unitAxis.X*inverseCos
	m[0][1] = unitAxis.X*unitAxis.Y*inverseCos + sin*unitAxis.Z
	m[0][2] = unitAxis.X*unitAxis.Z*inverseCos - sin*unitAxis.Y
	m[0][3] = 0
	m[1][0] = unitAxis.Y*unitAxis.X*inverseCos - sin*unitAxis.Z
	m[1][1] = cos + unitAxis.Y*unitAxis.Y*inverseCos
	m[1][2] = unitAxis.Y*unitAxis.Z*inverseCos + sin*unitAxis.X
	m[1][3] = 0
	m[2][0] = unitAxis.Z*unitAxis.X*inverseCos + sin*unitAxis.Y
	m[2][1] = unitAxis.Z*unitAxis.Y*inverseCos - sin*unitAxis.X
	m[2][2] = cos + unitAxis.Z*unitAxis.Z*inverseCos
	m[2][3] = 0
	m[3][0] = 0
	m[3][1] = 0
	m[3][2] = 0
	m[3][3] = 1

	return m
}

// RotateAroundAxis transforms this matrix by applying a rotation around the provided axis by the
// provided angle in radians
//
// axis - A 3-element vector that is normal to the angle of rotation. It does not need to be normalized.
//
// angleRad - The amount to rotate in radians
func (m *Mat4x4[T]) RotateAroundAxis(axis *Vec3[T], angleRad float64) *Mat4x4[T] {
	var unitAxis Vec3[T]
	unitAxis.SetVec3(axis).Normalize()

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

	m30 := rotate00*m[3][0] + rotate10*m[3][1] + rotate20*m[3][2]
	m31 := rotate01*m[3][0] + rotate11*m[3][1] + rotate21*m[3][2]
	m32 := rotate02*m[3][0] + rotate12*m[3][1] + rotate22*m[3][2]

	m[0][0] = m00
	m[0][1] = m01
	m[0][2] = m02
	m[1][0] = m10
	m[1][1] = m11
	m[1][2] = m12
	m[2][0] = m20
	m[2][1] = m21
	m[2][2] = m22
	m[3][0] = m30
	m[3][1] = m31
	m[3][2] = m32

	return m
}

// RotateX applies a transformation to this matrix that rotates around the x axis by the specified amount
//
// angleRad - The angle to rotate around the x axis in radians
func (m *Mat4x4[T]) RotateX(angleRad float64) *Mat4x4[T] {
	cos := T(math.Cos(angleRad))
	sin := T(math.Sin(angleRad))

	m01 := cos*m[0][1] - sin*m[0][2]
	m11 := cos*m[1][1] - sin*m[1][2]
	m21 := cos*m[2][1] - sin*m[2][2]
	m31 := cos*m[3][1] - sin*m[3][2]

	m02 := sin*m[0][1] + cos*m[0][2]
	m12 := sin*m[1][1] + cos*m[1][2]
	m22 := sin*m[2][1] + cos*m[2][2]
	m32 := sin*m[3][1] + cos*m[3][2]

	m[0][1] = m01
	m[1][1] = m11
	m[2][1] = m21
	m[3][1] = m31
	m[0][2] = m02
	m[1][2] = m12
	m[2][2] = m22
	m[3][2] = m32

	return m
}

// RotateY applies a transformation to this matrix that rotates around the y axis by the specified amount
//
// angleRad - The angle to rotate around the y axis in radians
func (m *Mat4x4[T]) RotateY(angleRad float64) *Mat4x4[T] {
	cos := T(math.Cos(angleRad))
	sin := T(math.Sin(angleRad))

	m00 := cos*m[0][0] + sin*m[0][2]
	m10 := cos*m[1][0] + sin*m[1][2]
	m20 := cos*m[2][0] + sin*m[2][2]
	m30 := cos*m[3][0] + sin*m[3][2]

	m02 := -sin*m[0][0] + cos*m[0][2]
	m12 := -sin*m[1][0] + cos*m[1][2]
	m22 := -sin*m[2][0] + cos*m[2][2]
	m32 := -sin*m[3][0] + cos*m[3][2]

	m[0][0] = m00
	m[1][0] = m10
	m[2][0] = m20
	m[3][0] = m30
	m[0][2] = m02
	m[1][2] = m12
	m[2][2] = m22
	m[3][2] = m32

	return m
}

// RotateZ applies a transformation to this matrix that rotates around the z axis by the specified amount
//
// angleRad - The angle to rotate around the z axis in radians
func (m *Mat4x4[T]) RotateZ(angleRad float64) *Mat4x4[T] {
	cos := T(math.Cos(angleRad))
	sin := T(math.Sin(angleRad))

	m00 := cos*m[0][0] - sin*m[0][1]
	m10 := cos*m[1][0] - sin*m[1][1]
	m20 := cos*m[2][0] - sin*m[2][1]
	m30 := cos*m[3][0] - sin*m[3][1]

	m01 := sin*m[0][0] + cos*m[0][1]
	m11 := sin*m[1][0] + cos*m[1][1]
	m21 := sin*m[2][0] + cos*m[2][1]
	m31 := sin*m[3][0] + cos*m[3][1]

	m[0][0] = m00
	m[1][0] = m10
	m[2][0] = m20
	m[3][0] = m30
	m[0][1] = m01
	m[1][1] = m11
	m[2][1] = m21
	m[3][1] = m31

	return m
}

// SetShearX overwrites this matrix with a shear matrix that shears along the X
// axis by the provided factor scalars
//
// y - y shear factor
//
// z - z shear factor
func (m *Mat4x4[T]) SetShearX(y, z T) *Mat4x4[T] {
	m[0][0] = 1
	m[1][0] = 0
	m[2][0] = 0
	m[3][0] = 0
	m[0][1] = y
	m[1][1] = 1
	m[2][1] = 0
	m[3][1] = 0
	m[0][2] = z
	m[1][2] = 0
	m[2][2] = 1
	m[3][2] = 0
	m[0][3] = 0
	m[1][3] = 0
	m[2][3] = 0
	m[3][3] = 1

	return m
}

// ShearX applies a shear transform to this matrix along the X axis by the provided
// factor scalars
//
// y - y shear factor
//
// z - z shear factor
func (m *Mat4x4[T]) ShearX(y, z T) *Mat4x4[T] {

	m01 := y*m[0][0] + m[0][1]
	m11 := y*m[1][0] + m[1][1]
	m21 := y*m[2][0] + m[2][1]
	m31 := y*m[3][0] + m[3][1]

	m02 := z*m[0][0] + m[0][2]
	m12 := z*m[1][0] + m[1][2]
	m22 := z*m[2][0] + m[2][2]
	m32 := z*m[3][0] + m[3][2]

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

// SetShearY overwrites this matrix with a shear matrix that shears along the Y
// axis by the provided factor scalars
//
// x - x shear factor
//
// z - z shear factor
func (m *Mat4x4[T]) SetShearY(x, z T) *Mat4x4[T] {
	m[0][0] = 1
	m[1][0] = x
	m[2][0] = 0
	m[3][0] = 0
	m[0][1] = 0
	m[1][1] = 1
	m[2][1] = 0
	m[3][1] = 0
	m[0][2] = 0
	m[1][2] = z
	m[2][2] = 1
	m[3][2] = 0
	m[0][3] = 0
	m[1][3] = 0
	m[2][3] = 0
	m[3][3] = 1

	return m
}

// ShearY applies a shear transform to this matrix along the Y axis by the provided
// factor scalars
//
// x - x shear factor
//
// z - z shear factor
func (m *Mat4x4[T]) ShearY(x, z T) *Mat4x4[T] {
	m00 := m[0][0] + x*m[0][1]
	m10 := m[1][0] + x*m[1][1]
	m20 := m[2][0] + x*m[2][1]
	m30 := m[3][0] + x*m[3][1]

	m02 := z*m[0][1] + m[0][2]
	m12 := z*m[1][1] + m[1][2]
	m22 := z*m[2][1] + m[2][2]
	m32 := z*m[3][1] + m[3][2]

	m[0][0] = m00
	m[0][2] = m02
	m[1][0] = m10
	m[1][2] = m12
	m[2][0] = m20
	m[2][2] = m22
	m[3][0] = m30
	m[3][2] = m32

	return m
}

// SetShearZ overwrites this matrix with a shear matrix that shears along the Z
// axis by the provided factor scalars
//
// x - x shear factor
//
// y - y shear factor
func (m *Mat4x4[T]) SetShearZ(x, y T) *Mat4x4[T] {
	m[0][0] = 1
	m[1][0] = 0
	m[2][0] = x
	m[3][0] = 0
	m[0][1] = 0
	m[1][1] = 1
	m[2][1] = y
	m[3][1] = 0
	m[0][2] = 0
	m[1][2] = 0
	m[2][2] = 1
	m[3][2] = 0
	m[0][3] = 0
	m[1][3] = 0
	m[2][3] = 0
	m[3][3] = 1

	return m
}

// ShearZ applies a shear transform to this matrix along the Z axis by the provided
// factor scalars
//
// x - x shear factor
//
// y - y shear factor
func (m *Mat4x4[T]) ShearZ(x, y T) *Mat4x4[T] {

	m00 := m[0][0] + x*m[0][2]
	m10 := m[1][0] + x*m[1][2]
	m20 := m[2][0] + x*m[2][2]
	m30 := m[3][0] + x*m[3][2]

	m01 := m[0][1] + y*m[0][2]
	m11 := m[1][1] + y*m[1][2]
	m21 := m[2][1] + y*m[2][2]
	m31 := m[3][1] + y*m[3][2]

	m[0][0] = m00
	m[0][1] = m01
	m[1][0] = m10
	m[1][1] = m11
	m[2][0] = m20
	m[2][1] = m21
	m[3][0] = m30
	m[3][1] = m31

	return m
}

// InterpolateMat4x4 updates this transform matrix by interpolating the rotation and translation
// with those of another transform matrix.
//
// otherMatrix - The "end" transform matrix in the interpolation
//
// delta - A value between 0 and 1 indicating how far to interpolate between this matrix and the other
func (m *Mat4x4[T]) InterpolateMat4x4(otherMatrix *Mat4x4[T], delta T) *Mat4x4[T] {
	var oldMatrix Mat4x4[T]
	oldMatrix.SetMat4x4(m)

	var thisRotation Mat4x4[T]
	thisRotation.SetAffineMat4x4(m)

	var transposedRotation Mat4x4[T]
	transposedRotation.SetMat4x4(&thisRotation).Transpose()

	var deltaRotation Mat4x4[T]
	deltaRotation.SetMultMat4x4(otherMatrix, &transposedRotation)

	var deltaAxis Vec3[T]
	var deltaAngle float64
	deltaRotation.GetAxisAngle(&deltaAxis, &deltaAngle)

	m.SetRotateAroundAxis(&deltaAxis, deltaAngle).MultMat4x4(&thisRotation)
	m[3][0] = oldMatrix[3][0] + delta*(otherMatrix[3][0]-oldMatrix[3][0])
	m[3][1] = oldMatrix[3][1] + delta*(otherMatrix[3][1]-oldMatrix[3][1])
	m[3][2] = oldMatrix[3][2] + delta*(otherMatrix[3][2]-oldMatrix[3][2])

	return m
}

// IsNormalized returns true if every row and every column in this matrix has a length of 1
//
// epsilon - The epsilon value to use in floating point comparisons. This much floating point
// drift is permitted before the method returns false. 0.0001 is a common epsilon value
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

// IsNull returns true if every column in this matrix has a length of 1
//
// epsilon - The epsilon value to use in floating point comparisons. This much floating point
// drift is permitted before the method returns false. 0.0001 is a common epsilon value
func (m *Mat4x4[T]) IsNull(epsilon T) bool {
	for i := 0; i < 4; i++ {
		column := abs[T](m[i][0]*m[i][0] + m[i][1]*m[i][1] + m[i][2]*m[i][2] + m[i][3]*m[i][3])
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
