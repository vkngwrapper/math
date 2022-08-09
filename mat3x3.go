package math

type Mat3Row[T FloatingPoint] [3]T
type Mat3x3[T FloatingPoint] [3]Mat3Row[T]

func (m *Mat3x3[T]) SetIdentity() *Mat3x3[T] {
	m[0][0] = 1
	m[0][1] = 0
	m[0][2] = 0
	m[1][0] = 0
	m[1][1] = 1
	m[1][2] = 0
	m[2][0] = 0
	m[2][1] = 0
	m[2][2] = 1

	return m
}

func (m *Mat3x3[T]) SetDiagonalScalar(s T) *Mat3x3[T] {
	m[0][0] = s
	m[0][1] = 0
	m[0][2] = 0
	m[1][0] = 0
	m[1][1] = s
	m[1][2] = 0
	m[2][0] = 0
	m[2][1] = 0
	m[2][2] = s

	return m
}

func (m *Mat3x3[T]) SetDiagonalVector(v *Vec3[T]) *Mat3x3[T] {
	m[0][0] = v.X
	m[0][1] = 0
	m[0][2] = 0
	m[1][0] = 0
	m[1][1] = v.Y
	m[1][2] = 0
	m[2][0] = 0
	m[2][1] = 0
	m[2][2] = v.Z

	return m
}

func (m *Mat3x3[T]) SetMat3x3(other *Mat3x3[T]) *Mat3x3[T] {
	m[0][0] = other[0][0]
	m[0][1] = other[0][1]
	m[0][2] = other[0][2]
	m[1][0] = other[1][0]
	m[1][1] = other[1][1]
	m[1][2] = other[1][2]
	m[2][0] = other[2][0]
	m[2][1] = other[2][1]
	m[2][2] = other[2][2]

	return m
}

func (m *Mat3x3[T]) SetMat4x4(other *Mat4x4[T]) *Mat3x3[T] {
	m[0][0] = other[0][0]
	m[0][1] = other[0][1]
	m[0][2] = other[0][2]
	m[1][0] = other[1][0]
	m[1][1] = other[1][1]
	m[1][2] = other[1][2]
	m[2][0] = other[2][0]
	m[2][1] = other[2][1]
	m[2][2] = other[2][2]

	return m
}

func (m *Mat3x3[T]) SetColMajor(c1, c2, c3 *Vec3[T]) *Mat3x3[T] {
	m[0][0] = c1.X
	m[0][1] = c1.Y
	m[0][2] = c1.Z
	m[1][0] = c2.X
	m[1][1] = c2.Y
	m[1][2] = c2.Z
	m[2][0] = c3.X
	m[2][1] = c3.Y
	m[2][2] = c3.Z

	return m
}

func (m *Mat3x3[T]) SetRowMajor(r1, r2, r3 *Vec3[T]) *Mat3x3[T] {
	m[0][0] = r1.X
	m[1][0] = r1.Y
	m[2][0] = r2.Z
	m[0][1] = r2.X
	m[1][1] = r2.Y
	m[2][1] = r2.Z
	m[0][2] = r3.X
	m[1][2] = r3.Y
	m[2][2] = r3.Z

	return m
}

func (m *Mat3x3[T]) SetMultMatrix(lhs, rhs *Mat3x3[T]) *Mat3x3[T] {
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

	return m
}

func (m *Mat3x3[T]) SetMatrixCrossProduct(vec *Vec3[T]) *Mat3x3[T] {
	m[0][0] = 0
	m[0][1] = vec.Z
	m[0][2] = -vec.Y
	m[1][0] = -vec.Z
	m[1][1] = 0
	m[1][2] = vec.X
	m[2][0] = vec.Y
	m[2][1] = -vec.X
	m[2][2] = 0

	return m
}

func (m *Mat3x3[T]) Transpose() *Mat3x3[T] {
	m[1][0], m[0][1] = m[0][1], m[1][0]
	m[2][0], m[0][2] = m[0][2], m[2][0]
	m[2][1], m[1][2] = m[1][2], m[2][1]

	return m
}

func (m *Mat3x3[T]) Inverse() *Mat3x3[T] {
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

	return m
}

func (m *Mat3x3[T]) Orthonormalize() *Mat3x3[T] {
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

	return m
}

func (m *Mat3x3[T]) MultMatrix(other *Mat3x3[T]) *Mat3x3[T] {
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

	return m
}

func (m *Mat3x3[T]) Proj2D(normal *Vec3[T]) *Mat3x3[T] {
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

	return m
}

func (m *Mat3x3[T]) ShearX(y T) *Mat3x3[T] {
	m[0][0] = m[0][0] + m[0][1]*y
	m[1][0] = m[1][0] + m[1][1]*y
	m[2][0] = m[2][0] + m[2][1]*y

	return m
}

func (m *Mat3x3[T]) ShearY(x T) *Mat3x3[T] {
	m[0][1] = m[0][0]*x + m[0][1]
	m[1][1] = m[1][0]*x + m[1][1]
	m[2][1] = m[2][0]*x + m[2][1]

	return m
}

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

func (m *Mat3x3[T]) IsNull(epsilon T) bool {
	for i := 0; i < 3; i++ {
		column := abs[T](m[i][0]*m[i][0] + m[i][1]*m[i][1] + m[i][2]*m[i][2])
		if column > epsilon {
			return false
		}
	}

	return true
}

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
