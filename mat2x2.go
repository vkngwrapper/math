package math

type Mat2Row[T FloatingPoint] [2]T

// Mat2x2 is a two-dimensional array of floating point-compatible values that can be used
// for 2x2 matrix arithmetic
type Mat2x2[T FloatingPoint] [2]Mat2Row[T]

// SetIdentity overwrites the current contents of the matrix with the identity matrix
func (m *Mat2x2[T]) SetIdentity() {
	m[0][0] = 1
	m[0][1] = 0
	m[1][0] = 0
	m[1][1] = 1
}

// SetDiagonalScalar overwrites the current contents of the matrix with an identity matrix,
// multiplied by the provided scalar value
//
// s - The scalar value to multiply the matrix entries by
func (m *Mat2x2[T]) SetDiagonalScalar(s T) {
	m[0][0] = s
	m[0][1] = 0
	m[1][0] = 0
	m[1][1] = s
}

// SetDiagonalVector overwrites the current contents of the matrix with an identity matrix,
// in which the 1 values are replaced with the contents of a provided vector, with each vector
// element corresponding to a column in the matrix
//
// v - The vector containing the diagonal entries that will be populated into the matrix
func (m *Mat2x2[T]) SetDiagonalVector(v *Vec2[T]) {
	m[0][0] = v.X
	m[0][1] = 0
	m[1][0] = 0
	m[1][1] = v.Y
}

// SetMat2x2 overwrites the current contents of the matrix with the contents of the provided
// matrix
//
// other - The matrix to initialize from
func (m *Mat2x2[T]) SetMat2x2(other *Mat2x2[T]) {
	m[0][0] = other[0][0]
	m[0][1] = other[0][1]
	m[1][0] = other[1][0]
	m[1][1] = other[1][1]
}

// SetMat3x3 overwrites the current contents of the matrix with the upper left
// 2x2 area of a 3x3 matrix
//
// other - The 3x3 matrix to initialize from
func (m *Mat2x2[T]) SetMat3x3(other *Mat3x3[T]) {
	m[0][0] = other[0][0]
	m[0][1] = other[0][1]
	m[1][0] = other[1][0]
	m[1][1] = other[1][1]
}

// SetMat4x4 overwrites the current contents of the matrix with the upper left
// 2x2 area of a 4x4 matrix
//
// other - The 4x4 matrix to initialize from
func (m *Mat2x2[T]) SetMat4x4(other *Mat4x4[T]) {
	m[0][0] = other[0][0]
	m[0][1] = other[0][1]
	m[1][0] = other[1][0]
	m[1][1] = other[1][1]
}

// SetColMajor overwrites the current contents of the matrix with 2 matrix columns
// passed in as 2-element vectors
//
// c1 - The first column of the matrix
// c2 - The second column of the matrix
func (m *Mat2x2[T]) SetColMajor(c1, c2 *Vec2[T]) {
	m[0][0] = c1.X
	m[0][1] = c1.Y
	m[1][0] = c2.X
	m[1][1] = c2.Y
}

// SetRowMajor overwrites the current contents of the matrix with 2 matrix rows
// passed in as 2-element vectors
//
// r1 - The first row of the matrix
// r2 - The second row of the matrix
func (m *Mat2x2[T]) SetRowMajor(r1, r2 *Vec2[T]) {
	m[0][0] = r1.X
	m[1][0] = r1.Y
	m[0][1] = r2.X
	m[1][1] = r2.Y
}

// SetMultMat2x2 multiplies two 2x2 matrices together and overwrites the current contents
// of this matrix with the results.
//
// lhs - The left operand of the multiplication operation
// rhs - The right operand of the multiplication operation
func (m *Mat2x2[T]) SetMultMat2x2(lhs, rhs *Mat2x2[T]) {
	m00 := lhs[0][0]*rhs[0][0] + lhs[1][0]*rhs[0][1]
	m10 := lhs[0][0]*rhs[1][0] + lhs[1][0]*rhs[1][1]

	m01 := lhs[0][1]*rhs[0][0] + lhs[1][1]*rhs[0][1]
	m11 := lhs[0][1]*rhs[1][0] + lhs[1][1]*rhs[1][1]

	m[0][0] = m00
	m[0][1] = m01
	m[1][0] = m10
	m[1][1] = m11
}

// MultMat2x2 multiplies this matrix against the provided matrix and updates this matrix
// with the results
//
// other - The right operand of the multiplication operation
func (m *Mat2x2[T]) MultMat2x2(other *Mat2x2[T]) {
	m00 := m[0][0]*other[0][0] + m[1][0]*other[0][1]
	m10 := m[0][0]*other[1][0] + m[1][0]*other[1][1]

	m01 := m[0][1]*other[0][0] + m[1][1]*other[0][1]
	m11 := m[0][1]*other[1][0] + m[1][1]*other[1][1]

	m[0][0] = m00
	m[0][1] = m01
	m[1][0] = m10
	m[1][1] = m11
}

// IsNormalized returns true if every row and every column in this matrix has a length of 1
//
// epsilon - The epsilon value to use in floating point comparisons. This much floating point
// drift is permitted before the method returns false. 0.0001 is a common epsilon value
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

// IsNull returns true if every column in this matrix has a length of 1
//
// epsilon - The epsilon value to use in floating point comparisons. This much floating point
// drift is permitted before the method returns false. 0.0001 is a common epsilon value
func (m *Mat2x2[T]) IsNull(epsilon T) bool {
	for i := 0; i < 2; i++ {
		column := abs[T](m[i][0]*m[i][0] + m[i][1]*m[i][1])
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
