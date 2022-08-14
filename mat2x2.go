package math

type Mat2Row [2]float32
type Mat2x2 [2]Mat2Row

func (m *Mat2x2) SetIdentity() *Mat2x2 {
	m[0][0] = 1
	m[0][1] = 0
	m[1][0] = 0
	m[1][1] = 1

	return m
}

func (m *Mat2x2) SetDiagonalScalar(s float32) *Mat2x2 {
	m[0][0] = s
	m[0][1] = 0
	m[1][0] = 0
	m[1][1] = s

	return m
}

func (m *Mat2x2) SetDiagonalVector(v *Vec2) *Mat2x2 {
	m[0][0] = v.X
	m[0][1] = 0
	m[1][0] = 0
	m[1][1] = v.Y

	return m
}

func (m *Mat2x2) SetMat2x2(other *Mat2x2) *Mat2x2 {
	m[0][0] = other[0][0]
	m[0][1] = other[0][1]
	m[1][0] = other[1][0]
	m[1][1] = other[1][1]

	return m
}

func (m *Mat2x2) SetMat3x3(other *Mat3x3) *Mat2x2 {
	m[0][0] = other[0][0]
	m[0][1] = other[0][1]
	m[1][0] = other[1][0]
	m[1][1] = other[1][1]

	return m
}

func (m *Mat2x2) SetMat4x4(other *Mat4x4) *Mat2x2 {
	m[0][0] = other[0][0]
	m[0][1] = other[0][1]
	m[1][0] = other[1][0]
	m[1][1] = other[1][1]

	return m
}

func (m *Mat2x2) SetColMajor(c1, c2 *Vec2) *Mat2x2 {
	m[0][0] = c1.X
	m[0][1] = c1.Y
	m[1][0] = c2.X
	m[1][1] = c2.Y

	return m
}

func (m *Mat2x2) SetRowMajor(r1, r2 *Vec2) *Mat2x2 {
	m[0][0] = r1.X
	m[1][0] = r1.Y
	m[0][1] = r2.X
	m[1][1] = r2.Y

	return m
}

func (m *Mat2x2) SetMultMatrix2x2(lhs, rhs *Mat2x2) *Mat2x2 {
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

func (m *Mat2x2) MultMatrix2x2(other *Mat2x2) *Mat2x2 {
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

func (m *Mat2x2) IsNormalized(epsilon float32) bool {
	for i := 0; i < 2; i++ {
		column := abs(m[i][0]*m[i][0] + m[i][1]*m[i][1])
		if column-1 > 2*epsilon {
			return false
		}
	}

	for i := 0; i < 2; i++ {
		row := abs(m[0][i]*m[0][i] + m[1][i]*m[1][i])
		for row-1 > 2*epsilon {
			return false
		}
	}

	return true
}

func (m *Mat2x2) IsNull(epsilon float32) bool {
	for i := 0; i < 2; i++ {
		column := abs(m[i][0]*m[i][0] + m[i][1]*m[i][1])
		if column > epsilon {
			return false
		}
	}

	return true
}

func (m *Mat2x2) Equal(other *Mat2x2, epsilon float32) bool {
	for i := 0; i < 2; i++ {
		for j := 0; j < 2; j++ {
			if abs(m[i][j]-other[i][j]) > epsilon {
				return false
			}
		}
	}

	return true
}
