package math

type Mat2Row[T FloatingPoint] [2]T
type Mat2x2[T FloatingPoint] [2]Mat2Row[T]

func (m *Mat2x2[T]) SetIdentity() *Mat2x2[T] {
	m[0][0] = 1
	m[0][1] = 0
	m[1][0] = 0
	m[1][1] = 1

	return m
}

func (m *Mat2x2[T]) SetDiagonalScalar(s T) *Mat2x2[T] {
	m[0][0] = s
	m[0][1] = 0
	m[1][0] = 0
	m[1][1] = s

	return m
}

func (m *Mat2x2[T]) SetDiagonalVector(v *Vec2[T]) *Mat2x2[T] {
	m[0][0] = v.X
	m[0][1] = 0
	m[1][0] = 0
	m[1][1] = v.Y

	return m
}

func (m *Mat2x2[T]) SetMat2x2(other *Mat2x2[T]) *Mat2x2[T] {
	m[0][0] = other[0][0]
	m[0][1] = other[0][1]
	m[1][0] = other[1][0]
	m[1][1] = other[1][1]

	return m
}

func (m *Mat2x2[T]) SetMat3x3(other *Mat3x3[T]) *Mat2x2[T] {
	m[0][0] = other[0][0]
	m[0][1] = other[0][1]
	m[1][0] = other[1][0]
	m[1][1] = other[1][1]

	return m
}

func (m *Mat2x2[T]) SetMat4x4(other *Mat4x4[T]) *Mat2x2[T] {
	m[0][0] = other[0][0]
	m[0][1] = other[0][1]
	m[1][0] = other[1][0]
	m[1][1] = other[1][1]

	return m
}

func (m *Mat2x2[T]) SetColMajor(c1, c2 *Vec2[T]) *Mat2x2[T] {
	m[0][0] = c1.X
	m[0][1] = c1.Y
	m[1][0] = c2.X
	m[1][1] = c2.Y

	return m
}

func (m *Mat2x2[T]) SetRowMajor(r1, r2 *Vec2[T]) *Mat2x2[T] {
	m[0][0] = r1.X
	m[1][0] = r1.Y
	m[0][1] = r2.X
	m[1][1] = r2.Y

	return m
}

func (m *Mat2x2[T]) SetMultMatrix2x2(lhs, rhs *Mat2x2[T]) *Mat2x2[T] {
	m00 := lhs[0][0]*rhs[0][0] + lhs[1][0]*rhs[0][1]
	m10 := lhs[0][0]*rhs[1][0] + lhs[1][0]*rhs[1][1]

	m01 := lhs[0][1]*rhs[0][0] + lhs[1][1]*rhs[0][1]
	m11 := lhs[0][1]*rhs[1][0] + lhs[1][1]*rhs[1][1]

	m[0][0] = m00
	m[0][1] = m01
	m[1][0] = m10
	m[1][1] = m11

	return m
}

func (m *Mat2x2[T]) MultMatrix2x2(other *Mat2x2[T]) *Mat2x2[T] {
	m00 := m[0][0]*other[0][0] + m[1][0]*other[0][1]
	m10 := m[0][0]*other[1][0] + m[1][0]*other[1][1]

	m01 := m[0][1]*other[0][0] + m[1][1]*other[0][1]
	m11 := m[0][1]*other[1][0] + m[1][1]*other[1][1]

	m[0][0] = m00
	m[0][1] = m01
	m[1][0] = m10
	m[1][1] = m11

	return m
}

func (m *Mat2x2[T]) IsNormalized(epsilon T) bool {
	for i := 0; i < 2; i++ {
		column := abs[T](m[i][0]*m[i][0] + m[i][1]*m[i][1])
		if column-T(1) > T(2)*epsilon {
			return false
		}
	}

	for i := 0; i < 2; i++ {
		row := abs[T](m[0][i]*m[0][i] + m[1][i]*m[1][i])
		for row-T(1) > T(2)*epsilon {
			return false
		}
	}

	return true
}

func (m *Mat2x2[T]) IsNull(epsilon T) bool {
	for i := 0; i < 2; i++ {
		column := abs[T](m[i][0]*m[i][0] + m[i][1]*m[i][1])
		if column > epsilon {
			return false
		}
	}

	return true
}

func (m *Mat2x2[T]) Equal(other *Mat2x2[T], epsilon T) bool {
	for i := 0; i < 2; i++ {
		for j := 0; j < 2; j++ {
			if abs[T](m[i][j]-other[i][j]) > epsilon {
				return false
			}
		}
	}

	return true
}
