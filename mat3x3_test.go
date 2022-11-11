package math

//func TestMat3x3_SetMultMatrix3x3(t *testing.T) {
//	var mat1 Mat3x3[float32]
//	mat1.SetIdentity()
//
//	mat1.SetRotationY(math.Pi / 2.0)
//
//	var mat2 Mat4x4[float32]
//	mat2.SetRotationY(math.Pi / 2.0)
//
//	var expected Mat4x4[float32]
//	expected.SetRotationY(math.Pi)
//
//	var result Mat4x4[float32]
//	result.SetMultMatrix4x4(&mat1, &mat2)
//
//	require.True(t, result.Equal(&expected, 0.0001))
//	require.False(t, mat1.Equal(&expected, 0.0001))
//	require.False(t, mat2.Equal(&expected, 0.0001))
//}
