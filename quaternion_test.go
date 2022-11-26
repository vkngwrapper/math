package math

import (
	"github.com/stretchr/testify/require"
	"math"
	"testing"
)

func TestQuaternion_RotateY_Clockwise(t *testing.T) {
	var quat Quaternion[float32]
	quat.SetIdentity()
	quat.RotateY(math.Pi / 2.0)

	vec := Vec3[float32]{1, 0, 0}
	vec.RotateWithQuaternion(&quat)

	require.InDelta(t, 0, vec.X, 0.0001)
	require.InDelta(t, 0, vec.Y, 0.0001)
	require.InDelta(t, -1, vec.Z, 0.0001)
}

func TestQuaternion_RotateY_CounterClockwise(t *testing.T) {
	var quat Quaternion[float32]
	quat.SetIdentity()
	quat.RotateY(-math.Pi / 2.0)

	vec := Vec3[float32]{1, 0, 0}
	vec.RotateWithQuaternion(&quat)

	require.InDelta(t, 0, vec.X, 0.0001)
	require.InDelta(t, 0, vec.Y, 0.0001)
	require.InDelta(t, 1, vec.Z, 0.0001)
}

func TestQuaternion_RotateX_Clockwise(t *testing.T) {
	var quat Quaternion[float32]
	quat.SetIdentity()
	quat.RotateX(math.Pi / 2.0)

	vec := Vec3[float32]{0, 0, -1}
	vec.RotateWithQuaternion(&quat)

	require.InDelta(t, 0, vec.X, 0.0001)
	require.InDelta(t, 1, vec.Y, 0.0001)
	require.InDelta(t, 0, vec.Z, 0.0001)
}

func TestQuaternion_RotateX_CounterClockwise(t *testing.T) {
	var quat Quaternion[float32]
	quat.SetIdentity()
	quat.RotateX(-math.Pi / 2.0)

	vec := Vec3[float32]{0, 0, -1}
	vec.RotateWithQuaternion(&quat)

	require.InDelta(t, 0, vec.X, 0.0001)
	require.InDelta(t, -1, vec.Y, 0.0001)
	require.InDelta(t, 0, vec.Z, 0.0001)
}

func TestQuaternion_RotateZ_Clockwise(t *testing.T) {
	var quat Quaternion[float32]
	quat.SetIdentity()
	quat.RotateZ(math.Pi / 2.0)

	vec := Vec3[float32]{1, 0, 0}
	vec.RotateWithQuaternion(&quat)

	require.InDelta(t, 0, vec.X, 0.0001)
	require.InDelta(t, 1, vec.Y, 0.0001)
	require.InDelta(t, 0, vec.Z, 0.0001)
}

func TestQuaternion_RotateZ_CounterClockwise(t *testing.T) {
	var quat Quaternion[float32]
	quat.SetIdentity()
	quat.RotateZ(-math.Pi / 2.0)

	vec := Vec3[float32]{1, 0, 0}
	vec.RotateWithQuaternion(&quat)

	require.InDelta(t, 0, vec.X, 0.0001)
	require.InDelta(t, -1, vec.Y, 0.0001)
	require.InDelta(t, 0, vec.Z, 0.0001)
}

func TestQuaternion_Rotate(t *testing.T) {
	var quat Quaternion[float32]
	quat.SetIdentity()
	quat.Rotate(&Vec3[float32]{1, 1, 0}, math.Pi)

	vec := Vec3[float32]{-1, 1, 0}
	vec.RotateWithQuaternion(&quat)

	require.InDelta(t, 1, vec.X, 0.0001)
	require.InDelta(t, -1, vec.Y, 0.0001)
	require.InDelta(t, 0, vec.Z, 0.0001)
}

func TestQuaternion_EulerAngles(t *testing.T) {
	var quat Quaternion[float32]
	quat.SetIdentity().RotateX(math.Pi / 2.0)
	require.InDelta(t, math.Pi/2.0, quat.Pitch(), 0.0001)

	quat.SetIdentity().RotateY(math.Pi / 3.0)
	require.InDelta(t, math.Pi/3.0, quat.Yaw(), 0.0001)

	quat.SetIdentity().RotateZ(math.Pi)
	require.InDelta(t, math.Pi, quat.Roll(), 0.0001)

	quat.SetIdentity().RotateX(math.Pi / 4.0).RotateZ(math.Pi)
	yaw, pitch, roll := quat.EulerAngles()
	require.InDelta(t, -math.Pi/4.0, pitch, 0.0001)
	require.InDelta(t, 0, yaw, 0.0001)
	require.InDelta(t, math.Pi, roll, 0.0001)
}

func TestQuaternion_SetLerp(t *testing.T) {
	var quat1 Quaternion[float32]
	quat1.SetRotateY(math.Pi)

	var quat2 Quaternion[float32]
	quat2.SetQuaternion(&quat1).RotateZ(math.Pi)

	start := Vec3[float32]{X: -1, Y: 0, Z: 0}

	var afterQuat1 Vec3[float32]
	afterQuat1.SetVec3(&start).RotateWithQuaternion(&quat1)
	require.InDelta(t, 1, afterQuat1.X, 0.0001)
	require.InDelta(t, 0, afterQuat1.Y, 0.0001)
	require.InDelta(t, 0, afterQuat1.Z, 0.0001)

	var afterQuat2 Vec3[float32]
	afterQuat2.SetVec3(&start).RotateWithQuaternion(&quat2)
	require.InDelta(t, -1, afterQuat2.X, 0.0001)
	require.InDelta(t, 0, afterQuat2.Y, 0.0001)
	require.InDelta(t, 0, afterQuat2.Z, 0.0001)

	var halfWayQuat Quaternion[float32]
	halfWayQuat.SetLerp(&quat1, &quat2, 0.5)

	var halfWay Vec3[float32]
	halfWay.SetVec3(&start).RotateWithQuaternion(&halfWayQuat)
	require.InDelta(t, 0, halfWay.X, 0.0001)
	require.InDelta(t, -1, halfWay.Y, 0.0001)
	require.InDelta(t, 0, halfWay.Z, 0.0001)
}

func TestQuaternion_Lerp(t *testing.T) {
	var quat1 Quaternion[float32]
	quat1.SetRotateY(math.Pi)

	var quat2 Quaternion[float32]
	quat2.SetQuaternion(&quat1).RotateZ(math.Pi)

	start := Vec3[float32]{X: -1, Y: 0, Z: 0}

	var afterQuat1 Vec3[float32]
	afterQuat1.SetVec3(&start).RotateWithQuaternion(&quat1)
	require.InDelta(t, 1, afterQuat1.X, 0.0001)
	require.InDelta(t, 0, afterQuat1.Y, 0.0001)
	require.InDelta(t, 0, afterQuat1.Z, 0.0001)

	var afterQuat2 Vec3[float32]
	afterQuat2.SetVec3(&start).RotateWithQuaternion(&quat2)
	require.InDelta(t, -1, afterQuat2.X, 0.0001)
	require.InDelta(t, 0, afterQuat2.Y, 0.0001)
	require.InDelta(t, 0, afterQuat2.Z, 0.0001)

	var halfWayQuat Quaternion[float32]
	halfWayQuat.SetQuaternion(&quat1).Lerp(&quat2, 0.5)

	var halfWay Vec3[float32]
	halfWay.SetVec3(&start).RotateWithQuaternion(&halfWayQuat)
	require.InDelta(t, 0, halfWay.X, 0.0001)
	require.InDelta(t, -1, halfWay.Y, 0.0001)
	require.InDelta(t, 0, halfWay.Z, 0.0001)
}

func TestQuaternion_SetVectorRotation_Corner(t *testing.T) {
	origin := Vec3[float32]{X: 1, Y: 0, Z: 0}
	dest := Vec3[float32]{X: 0, Y: 0, Z: 1}

	var quat Quaternion[float32]
	quat.SetVectorRotation(&origin, &dest)

	var replayDest Vec3[float32]
	replayDest.SetRotateWithQuaternion(&origin, &quat)
	require.InDelta(t, 0, replayDest.X, 0.0001)
	require.InDelta(t, 0, replayDest.Y, 0.0001)
	require.InDelta(t, 1, replayDest.Z, 0.0001)
}

func TestQuaternion_SetVectorRotation_Line(t *testing.T) {
	origin := Vec3[float32]{X: 1, Y: 0, Z: 0}
	dest := Vec3[float32]{X: -1, Y: 0, Z: 0}

	var quat Quaternion[float32]
	quat.SetVectorRotation(&origin, &dest)

	var replayDest Vec3[float32]
	replayDest.SetRotateWithQuaternion(&origin, &quat)
	require.InDelta(t, -1, replayDest.X, 0.0001)
	require.InDelta(t, 0, replayDest.Y, 0.0001)
	require.InDelta(t, 0, replayDest.Z, 0.0001)
}

func TestQuaternion_Inverse(t *testing.T) {
	origin := Vec3[float32]{X: 1, Y: 0, Z: 0}
	dest := Vec3[float32]{X: 0, Y: 0, Z: 1}

	var quat Quaternion[float32]
	quat.SetVectorRotation(&origin, &dest).Inverse()

	var replayOrigin Vec3[float32]
	replayOrigin.SetRotateWithQuaternion(&dest, &quat)
	require.InDelta(t, 1, replayOrigin.X, 0.0001)
	require.InDelta(t, 0, replayOrigin.Y, 0.0001)
	require.InDelta(t, 0, replayOrigin.Z, 0.0001)
}

func TestQuaternion_SetAxisAngle(t *testing.T) {
	axis := Vec3[float32]{X: 0, Y: 1, Z: 0}
	var quat Quaternion[float32]
	quat.SetAxisAngle(&axis, math.Pi)

	vec := Vec3[float32]{X: 1, Y: 0, Z: 0}
	vec.RotateWithQuaternion(&quat)
	require.InDelta(t, -1, vec.X, 0.0001)
	require.InDelta(t, 0, vec.Y, 0.0001)
	require.InDelta(t, 0, vec.Z, 0.0001)

	var outputAxis Vec3[float32]
	quat.GetAxis(&outputAxis)
	require.InDelta(t, 0, outputAxis.X, 0.0001)
	require.InDelta(t, 1, outputAxis.Y, 0.0001)
	require.InDelta(t, 0, outputAxis.Z, 0.0001)

	outputAngle := quat.Angle()
	require.InDelta(t, math.Pi, outputAngle, 0.0001)
}

func TestQuaternion_DotProduct(t *testing.T) {
	var quat1 Quaternion[float32]
	quat1.SetRotateX(math.Pi / 2.0).Normalize()

	var quat2 Quaternion[float32]
	quat2.SetRotateX(3.0 * math.Pi / 2.0).Normalize()

	cosDiff := quat1.DotProduct(&quat2)
	require.InDelta(t, math.Pi, math.Acos(float64(cosDiff))*2.0, 0.0001)
}

func TestQuaternion_SetMat3x3(t *testing.T) {
	var rotMat Mat3x3[float32]
	rotMat.SetRotationY(math.Pi)

	var quat Quaternion[float32]
	quat.SetMat3x3(&rotMat)

	var outMat Mat3x3[float32]
	outMat.SetQuaternion(&quat)

	require.True(t, rotMat.Equal(&outMat, 0.0001))

	matStart := Vec3[float32]{X: 1, Y: 0, Z: 0}
	matStart.Transform(&rotMat)

	quatStart := Vec3[float32]{X: 1, Y: 0, Z: 0}
	quatStart.RotateWithQuaternion(&quat)

	require.Equal(t, matStart, quatStart)
}

func TestQuaternion_SetMat4x4(t *testing.T) {
	var rotMat Mat4x4[float32]
	rotMat.SetRotationY(math.Pi)

	var quat Quaternion[float32]
	quat.SetMat4x4(&rotMat)

	var outMat Mat4x4[float32]
	outMat.SetQuaternion(&quat)

	require.True(t, rotMat.Equal(&outMat, 0.0001))

	matStart := Vec4[float32]{X: 1, Y: 0, Z: 0, W: 1}
	matStart.Transform(&rotMat)

	quatStart := Vec4[float32]{X: 1, Y: 0, Z: 0, W: 1}
	quatStart.RotateWithQuaternion(&quat)

	require.Equal(t, matStart, quatStart)
}

func TestQuaternion_MultQuaternion(t *testing.T) {
	var q1, q2, q3 Quaternion[float32]
	q1.SetRotateX(math.Pi / 4.0)
	q2.SetRotateX(math.Pi / 4.0)
	q3.SetRotateX(math.Pi / 2.0)

	q1.MultQuaternion(&q2)

	require.True(t, q1.Equal(&q3, 0.0001))
}

func TestQuaternion_MultQuaternion_TwoPart(t *testing.T) {
	var q1, q2, mult, expected Quaternion[float32]
	q1.SetRotateY(math.Pi / 2.0)
	q2.SetRotateX(math.Pi / 2.0)

	mult.SetMultQuaternion(&q1, &q2)

	expected.SetRotateY(math.Pi / 2.0).RotateX(math.Pi / 2.0)

	require.True(t, expected.Equal(&mult, 0.0001))
}

func TestQuaternion_SetEulerAngles(t *testing.T) {
	var eulers Quaternion[float32]
	eulers.SetEulerAngles(math.Pi/2.0, math.Pi/4.0, math.Pi)

	vec := Vec3[float32]{X: 1, Y: 0, Z: 0}
	vec.RotateWithQuaternion(&eulers)

	require.InDelta(t, 0, vec.X, 0.0001)
	require.InDelta(t, float32(math.Sqrt(2.0)/2.0), vec.Y, 0.0001)
	require.InDelta(t, float32(-math.Sqrt(2.0)/2.0), vec.Z, 0.0001)
}

func TestQuaternion_SetConjugate(t *testing.T) {
	var eulers Quaternion[float32]
	eulers.SetEulerAngles(math.Pi/2.0, math.Pi/4.0, math.Pi)

	var inverse Quaternion[float32]
	inverse.SetConjugate(&eulers)

	var identity Quaternion[float32]
	identity.SetMultQuaternion(&eulers, &inverse)

	require.InDelta(t, 0.0, identity.X, 0.0001)
	require.InDelta(t, 0.0, identity.Y, 0.0001)
	require.InDelta(t, 0.0, identity.Z, 0.0001)
	require.InDelta(t, 1.0, identity.W, 0.0001)
}

func TestQuaternion_Conjugate(t *testing.T) {
	var eulers Quaternion[float32]
	eulers.SetEulerAngles(math.Pi/2.0, math.Pi/4.0, math.Pi)

	var inverse Quaternion[float32]
	inverse.SetQuaternion(&eulers).Conjugate()

	var identity Quaternion[float32]
	identity.SetMultQuaternion(&eulers, &inverse)

	require.InDelta(t, 0.0, identity.X, 0.0001)
	require.InDelta(t, 0.0, identity.Y, 0.0001)
	require.InDelta(t, 0.0, identity.Z, 0.0001)
	require.InDelta(t, 1.0, identity.W, 0.0001)
}

func TestQuaternion_Log(t *testing.T) {
	var yRot Quaternion[float32]
	yRot.SetRotateY(math.Pi).Log()

	require.InDelta(t, 0.0, yRot.X, 0.0001)
	require.InDelta(t, math.Pi/2.0, yRot.Y, 0.0001)
	require.InDelta(t, 0, yRot.Z, 0.0001)
	require.InDelta(t, 0, yRot.W, 0.0001)
}

func TestQuaternion_Exp(t *testing.T) {
	var yRot Quaternion[float32]
	yRot.SetRotateY(math.Pi).Log().Exp()

	vec := Vec3[float32]{X: 1, Y: 0, Z: 0}
	vec.RotateWithQuaternion(&yRot)

	require.InDelta(t, -1, vec.X, 0.0001)
	require.InDelta(t, 0, vec.Y, 0.0001)
	require.InDelta(t, 0, vec.Z, 0.0001)
}

func TestQuaternion_Pow_Zero(t *testing.T) {
	var yRot Quaternion[float32]
	yRot.SetRotateY(math.Pi).Pow(0)

	require.InDelta(t, 0.0, yRot.X, 0.0001)
	require.InDelta(t, 0.0, yRot.Y, 0.0001)
	require.InDelta(t, 0.0, yRot.Z, 0.0001)
	require.InDelta(t, 1.0, yRot.W, 0.0001)
}

func TestQuaternion_Pow_One(t *testing.T) {
	var yRot Quaternion[float32]
	yRot.SetRotateY(math.Pi).Pow(1)

	vec := Vec3[float32]{X: 1, Y: 0, Z: 0}
	vec.RotateWithQuaternion(&yRot)

	require.InDelta(t, -1, vec.X, 0.0001)
	require.InDelta(t, 0, vec.Y, 0.0001)
	require.InDelta(t, 0, vec.Z, 0.0001)
}

func TestQuaternion_Pow_Two(t *testing.T) {
	var yRot Quaternion[float32]
	yRot.SetRotateY(math.Pi).Pow(2)

	vec := Vec3[float32]{X: 1, Y: 0, Z: 0}
	vec.RotateWithQuaternion(&yRot)

	require.InDelta(t, 1, vec.X, 0.0001)
	require.InDelta(t, 0, vec.Y, 0.0001)
	require.InDelta(t, 0, vec.Z, 0.0001)
}

func TestQuaternion_Slerp(t *testing.T) {
	var quat1 Quaternion[float32]
	quat1.SetRotateY(math.Pi)

	var quat2 Quaternion[float32]
	quat2.SetQuaternion(&quat1).RotateZ(math.Pi)

	start := Vec3[float32]{X: -1, Y: 0, Z: 0}

	var afterQuat1 Vec3[float32]
	afterQuat1.SetVec3(&start).RotateWithQuaternion(&quat1)
	require.InDelta(t, 1, afterQuat1.X, 0.0001)
	require.InDelta(t, 0, afterQuat1.Y, 0.0001)
	require.InDelta(t, 0, afterQuat1.Z, 0.0001)

	var afterQuat2 Vec3[float32]
	afterQuat2.SetVec3(&start).RotateWithQuaternion(&quat2)
	require.InDelta(t, -1, afterQuat2.X, 0.0001)
	require.InDelta(t, 0, afterQuat2.Y, 0.0001)
	require.InDelta(t, 0, afterQuat2.Z, 0.0001)

	var halfWayQuat Quaternion[float32]
	halfWayQuat.SetQuaternion(&quat1).Slerp(&quat2, 0.5)

	var halfWay Vec3[float32]
	halfWay.SetVec3(&start).RotateWithQuaternion(&halfWayQuat)
	require.InDelta(t, 0, halfWay.X, 0.0001)
	require.InDelta(t, -1, halfWay.Y, 0.0001)
	require.InDelta(t, 0, halfWay.Z, 0.0001)
}

func TestQuaternion_SetSlerp(t *testing.T) {
	var quat1 Quaternion[float32]
	quat1.SetRotateY(math.Pi)

	var quat2 Quaternion[float32]
	quat2.SetQuaternion(&quat1).RotateZ(math.Pi)

	start := Vec3[float32]{X: -1, Y: 0, Z: 0}

	var afterQuat1 Vec3[float32]
	afterQuat1.SetVec3(&start).RotateWithQuaternion(&quat1)
	require.InDelta(t, 1, afterQuat1.X, 0.0001)
	require.InDelta(t, 0, afterQuat1.Y, 0.0001)
	require.InDelta(t, 0, afterQuat1.Z, 0.0001)

	var afterQuat2 Vec3[float32]
	afterQuat2.SetVec3(&start).RotateWithQuaternion(&quat2)
	require.InDelta(t, -1, afterQuat2.X, 0.0001)
	require.InDelta(t, 0, afterQuat2.Y, 0.0001)
	require.InDelta(t, 0, afterQuat2.Z, 0.0001)

	var halfWayQuat Quaternion[float32]
	halfWayQuat.SetSlerp(&quat1, &quat2, 0.5)

	var halfWay Vec3[float32]
	halfWay.SetVec3(&start).RotateWithQuaternion(&halfWayQuat)
	require.InDelta(t, 0, halfWay.X, 0.0001)
	require.InDelta(t, -1, halfWay.Y, 0.0001)
	require.InDelta(t, 0, halfWay.Z, 0.0001)
}

func TestQuaternion_Mix(t *testing.T) {
	var quat1 Quaternion[float32]
	quat1.SetRotateY(math.Pi)

	var quat2 Quaternion[float32]
	quat2.SetQuaternion(&quat1).RotateZ(math.Pi)

	start := Vec3[float32]{X: -1, Y: 0, Z: 0}

	var afterQuat1 Vec3[float32]
	afterQuat1.SetVec3(&start).RotateWithQuaternion(&quat1)
	require.InDelta(t, 1, afterQuat1.X, 0.0001)
	require.InDelta(t, 0, afterQuat1.Y, 0.0001)
	require.InDelta(t, 0, afterQuat1.Z, 0.0001)

	var afterQuat2 Vec3[float32]
	afterQuat2.SetVec3(&start).RotateWithQuaternion(&quat2)
	require.InDelta(t, -1, afterQuat2.X, 0.0001)
	require.InDelta(t, 0, afterQuat2.Y, 0.0001)
	require.InDelta(t, 0, afterQuat2.Z, 0.0001)

	var halfWayQuat Quaternion[float32]
	halfWayQuat.SetQuaternion(&quat1).Mix(&quat2, 0.5)

	var halfWay Vec3[float32]
	halfWay.SetVec3(&start).RotateWithQuaternion(&halfWayQuat)
	require.InDelta(t, 0, halfWay.X, 0.0001)
	require.InDelta(t, -1, halfWay.Y, 0.0001)
	require.InDelta(t, 0, halfWay.Z, 0.0001)
}

func TestQuaternion_SetMix(t *testing.T) {
	var quat1 Quaternion[float32]
	quat1.SetRotateY(math.Pi)

	var quat2 Quaternion[float32]
	quat2.SetQuaternion(&quat1).RotateZ(math.Pi)

	start := Vec3[float32]{X: -1, Y: 0, Z: 0}

	var afterQuat1 Vec3[float32]
	afterQuat1.SetVec3(&start).RotateWithQuaternion(&quat1)
	require.InDelta(t, 1, afterQuat1.X, 0.0001)
	require.InDelta(t, 0, afterQuat1.Y, 0.0001)
	require.InDelta(t, 0, afterQuat1.Z, 0.0001)

	var afterQuat2 Vec3[float32]
	afterQuat2.SetVec3(&start).RotateWithQuaternion(&quat2)
	require.InDelta(t, -1, afterQuat2.X, 0.0001)
	require.InDelta(t, 0, afterQuat2.Y, 0.0001)
	require.InDelta(t, 0, afterQuat2.Z, 0.0001)

	var halfWayQuat Quaternion[float32]
	halfWayQuat.SetMix(&quat1, &quat2, 0.5)

	var halfWay Vec3[float32]
	halfWay.SetVec3(&start).RotateWithQuaternion(&halfWayQuat)
	require.InDelta(t, 0, halfWay.X, 0.0001)
	require.InDelta(t, -1, halfWay.Y, 0.0001)
	require.InDelta(t, 0, halfWay.Z, 0.0001)
}

func TestQuaternion_ShortMix(t *testing.T) {
	var quat1 Quaternion[float32]
	quat1.SetRotateY(math.Pi)

	var quat2 Quaternion[float32]
	quat2.SetQuaternion(&quat1).RotateZ(math.Pi)

	start := Vec3[float32]{X: -1, Y: 0, Z: 0}

	var afterQuat1 Vec3[float32]
	afterQuat1.SetVec3(&start).RotateWithQuaternion(&quat1)
	require.InDelta(t, 1, afterQuat1.X, 0.0001)
	require.InDelta(t, 0, afterQuat1.Y, 0.0001)
	require.InDelta(t, 0, afterQuat1.Z, 0.0001)

	var afterQuat2 Vec3[float32]
	afterQuat2.SetVec3(&start).RotateWithQuaternion(&quat2)
	require.InDelta(t, -1, afterQuat2.X, 0.0001)
	require.InDelta(t, 0, afterQuat2.Y, 0.0001)
	require.InDelta(t, 0, afterQuat2.Z, 0.0001)

	var halfWayQuat Quaternion[float32]
	halfWayQuat.SetQuaternion(&quat1).ShortMix(&quat2, 0.5)

	var halfWay Vec3[float32]
	halfWay.SetVec3(&start).RotateWithQuaternion(&halfWayQuat)
	require.InDelta(t, 0, halfWay.X, 0.0001)
	require.InDelta(t, -1, halfWay.Y, 0.0001)
	require.InDelta(t, 0, halfWay.Z, 0.0001)
}

func TestQuaternion_SetShortMix(t *testing.T) {
	var quat1 Quaternion[float32]
	quat1.SetRotateY(math.Pi)

	var quat2 Quaternion[float32]
	quat2.SetQuaternion(&quat1).RotateZ(math.Pi)

	start := Vec3[float32]{X: -1, Y: 0, Z: 0}

	var afterQuat1 Vec3[float32]
	afterQuat1.SetVec3(&start).RotateWithQuaternion(&quat1)
	require.InDelta(t, 1, afterQuat1.X, 0.0001)
	require.InDelta(t, 0, afterQuat1.Y, 0.0001)
	require.InDelta(t, 0, afterQuat1.Z, 0.0001)

	var afterQuat2 Vec3[float32]
	afterQuat2.SetVec3(&start).RotateWithQuaternion(&quat2)
	require.InDelta(t, -1, afterQuat2.X, 0.0001)
	require.InDelta(t, 0, afterQuat2.Y, 0.0001)
	require.InDelta(t, 0, afterQuat2.Z, 0.0001)

	var halfWayQuat Quaternion[float32]
	halfWayQuat.SetShortMix(&quat1, &quat2, 0.5)

	var halfWay Vec3[float32]
	halfWay.SetVec3(&start).RotateWithQuaternion(&halfWayQuat)
	require.InDelta(t, 0, halfWay.X, 0.0001)
	require.InDelta(t, -1, halfWay.Y, 0.0001)
	require.InDelta(t, 0, halfWay.Z, 0.0001)
}

func TestQuaternion_Squad(t *testing.T) {
	var key1, key2, key3, key4, control3, control4 Quaternion[float32]
	key1.SetIdentity()
	key2.SetRotateY(2 * math.Pi / 3.0)
	key3.SetRotateY(4 * math.Pi / 3.0)
	key4.SetIdentity()

	control3.SetIntermediate(&key2, &key3, &key4)
	control4.SetIntermediate(&key3, &key4, &key4)

	var halfWayQuat Quaternion[float32]
	halfWayQuat.SetSquad(&key3, &key4, &control3, &control4, 1)

	vec := Vec3[float32]{X: 1, Y: 0, Z: 0}
	vec.RotateWithQuaternion(&halfWayQuat)

	require.InDelta(t, 1, vec.X, 0.0001)
	require.InDelta(t, 0, vec.Y, 0.0001)
	require.InDelta(t, 0, vec.Z, 0.0001)
}
