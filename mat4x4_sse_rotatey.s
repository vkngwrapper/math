#include "textflag.h"

// func __m4x4RotateY(m *Mat4x4, sin, cos float32)
TEXT ·__m4x4RotateY(SB),NOSPLIT,$0-16
    MOVQ m+0(FP), AX

    CALL math·Sin(SB)
    // Load cos & sin into the lower float of X0 and X1 registers
    MOVSS sin+8(FP), X0
    MOVSS cos+12(FP), X1

    // Broadcast the sin/cos value through the entire X0/X1 registers
    SHUFPS $0x0, X0, X0 // sin
    SHUFPS $0x0, X1, X1 // cos

    // Load rows into X registers
    MOVUPS 0(AX), X3
    MOVUPS 32(AX), X4

    //m00 := m[0][0]*cos - m[2][0]*sin
    //m01 := m[0][1]*cos - m[2][1]*sin
    //m02 := m[0][2]*cos - m[2][2]*sin
    //m03 := m[0][3]*cos - m[2][3]*sin

    MOVUPS X3, X5 // x5 = row 0
    MULPS X1, X5 // x5 *= cos
    MOVUPS X4, X6 // x6 = row 2
    MULPS X0, X6 // x6 *= cos
    SUBPS X6, X5 // x5 -= x6

    //m20 := m[0][0]*sin + m[2][0]*cos
    //m21 := m[0][1]*sin + m[2][1]*cos
    //m22 := m[0][2]*sin + m[2][2]*cos
    //m23 := m[0][3]*sin + m[2][3]*cos
    MOVUPS X3, X6 // x6 = row 0
    MULPS X0, X6 // x6 *= sin
    MOVUPS X4, X7 // x7 = row 1
    MULPS X1, X7 // x7 *= cos
    ADDPS X7, X6 // x6 += x7

    MOVUPS X5, m+0(FP) // row 0 = x5
    MOVUPS X6, m+32(FP) // row 2 = x6
    RET
