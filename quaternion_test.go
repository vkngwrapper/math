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

	quat.SetIdentity().RotateX(math.Pi / 4.0).RotateZ(math.Pi).RotateY(math.Pi / 3.0)
	yaw, pitch, roll := quat.EulerAngles()
	require.InDelta(t, -math.Pi/4.0, pitch, 0.0001)
	require.InDelta(t, math.Pi/3.0, yaw, 0.0001)
	require.InDelta(t, 0, roll, 0.0001)
}
