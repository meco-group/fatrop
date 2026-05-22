//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//
// Small header-only collection of quaternion / SO(3) helpers shared by the
// Lie-group examples. Everything lives in fatrop::examples to keep things out
// of the global namespace.
//

#ifndef __fatrop_examples_so3_helpers_hpp__
#define __fatrop_examples_so3_helpers_hpp__

#include "fatrop/context/context.hpp"
#include <cmath>

namespace fatrop
{
namespace examples
{

// ---- 3-vector helpers ------------------------------------------------------

inline void cross3(const Scalar *a, const Scalar *b, Scalar *r)
{
    r[0] = a[1] * b[2] - a[2] * b[1];
    r[1] = a[2] * b[0] - a[0] * b[2];
    r[2] = a[0] * b[1] - a[1] * b[0];
}

// 3x3 (row-major) helpers ----------------------------------------------------

inline void mat3_identity(Scalar M[9])
{
    for (int i = 0; i < 9; ++i) M[i] = 0.;
    M[0] = M[4] = M[8] = 1.;
}

inline void mat3_mul(const Scalar A[9], const Scalar B[9], Scalar C[9])
{
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
            Scalar s = 0.;
            for (int k = 0; k < 3; ++k) s += A[i * 3 + k] * B[k * 3 + j];
            C[i * 3 + j] = s;
        }
}

// Skew-symmetric "hat" map [v]_x.
inline void hat3(const Scalar v[3], Scalar M[9])
{
    M[0] = 0.;    M[1] = -v[2]; M[2] = v[1];
    M[3] = v[2];  M[4] = 0.;    M[5] = -v[0];
    M[6] = -v[1]; M[7] = v[0];  M[8] = 0.;
}

// ---- Unit quaternion helpers ([w, x, y, z] / Hamilton convention) ----------

// Hamilton quaternion product: r = a (x) b.
inline void quat_mul(const Scalar *a, const Scalar *b, Scalar *r)
{
    const Scalar aw = a[0], ax = a[1], ay = a[2], az = a[3];
    const Scalar bw = b[0], bx = b[1], by = b[2], bz = b[3];
    r[0] = aw * bw - ax * bx - ay * by - az * bz;
    r[1] = aw * bx + ax * bw + ay * bz - az * by;
    r[2] = aw * by - ax * bz + ay * bw + az * bx;
    r[3] = aw * bz + ax * by - ay * bx + az * bw;
}

// In-place renormalization for numerical robustness.
inline void quat_normalize(Scalar *q)
{
    const Scalar n = std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
    if (n > 0.)
    {
        q[0] /= n;
        q[1] /= n;
        q[2] /= n;
        q[3] /= n;
    }
}

// Inverse of a unit quaternion (conjugate).
inline void quat_inv(const Scalar *q, Scalar *qi)
{
    qi[0] = q[0];
    qi[1] = -q[1];
    qi[2] = -q[2];
    qi[3] = -q[3];
}

// Exponential map from so(3) to the unit quaternions:
//     Exp_q(v) = cos(|v|/2) + sin(|v|/2) * v / |v|
inline void quat_exp(const Scalar v[3], Scalar q[4])
{
    const Scalar half = 0.5 * std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (half < 1e-12)
    {
        // Second-order Taylor; avoids 0/0 at the identity.
        q[0] = 1. - 0.5 * half * half;
        q[1] = 0.5 * v[0];
        q[2] = 0.5 * v[1];
        q[3] = 0.5 * v[2];
        quat_normalize(q);
        return;
    }
    const Scalar s = std::sin(half) / (2. * half);
    q[0] = std::cos(half);
    q[1] = s * v[0];
    q[2] = s * v[1];
    q[3] = s * v[2];
}

// Logarithm on the unit quaternions: returns the rotation vector v in R^3
// such that Exp_q(v) = q (sign chosen so |v| in [0, pi]).
inline void quat_log(const Scalar *q, Scalar *v)
{
    Scalar w = q[0], x = q[1], y = q[2], z = q[3];
    if (w < 0.)
    {
        w = -w;
        x = -x;
        y = -y;
        z = -z;
    }
    const Scalar vec_norm = std::sqrt(x * x + y * y + z * z);
    if (vec_norm < 1e-12)
    {
        v[0] = 2. * x;
        v[1] = 2. * y;
        v[2] = 2. * z;
        return;
    }
    const Scalar angle = 2. * std::atan2(vec_norm, w);
    const Scalar s = angle / vec_norm;
    v[0] = s * x;
    v[1] = s * y;
    v[2] = s * z;
}

// Rotate a 3-vector p by a unit quaternion q: p_rot = R(q) p.
inline void quat_rotate(const Scalar *q, const Scalar *p, Scalar *p_rot)
{
    const Scalar w = q[0], x = q[1], y = q[2], z = q[3];
    const Scalar vp = x * p[0] + y * p[1] + z * p[2];
    const Scalar vv = x * x + y * y + z * z;
    const Scalar c = 1.0 - 2.0 * vv; // = w^2 - |v|^2 for unit q
    const Scalar cx = y * p[2] - z * p[1];
    const Scalar cy = z * p[0] - x * p[2];
    const Scalar cz = x * p[1] - y * p[0];
    p_rot[0] = c * p[0] + 2.0 * vp * x + 2.0 * w * cx;
    p_rot[1] = c * p[1] + 2.0 * vp * y + 2.0 * w * cy;
    p_rot[2] = c * p[2] + 2.0 * vp * z + 2.0 * w * cz;
}

// Convert a unit quaternion to a row-major 3x3 rotation matrix.
inline void quat_to_rotmat(const Scalar *q, Scalar R[9])
{
    const Scalar w = q[0], x = q[1], y = q[2], z = q[3];
    R[0] = 1 - 2 * (y * y + z * z); R[1] = 2 * (x * y - z * w); R[2] = 2 * (x * z + y * w);
    R[3] = 2 * (x * y + z * w);     R[4] = 1 - 2 * (x * x + z * z); R[5] = 2 * (y * z - x * w);
    R[6] = 2 * (x * z - y * w);     R[7] = 2 * (y * z + x * w); R[8] = 1 - 2 * (x * x + y * y);
}

// SO(3) left Jacobian J_l(v) (closed-form, small-angle Taylor near v = 0).
inline void so3_left_jacobian(const Scalar v[3], Scalar J[9])
{
    const Scalar theta2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    mat3_identity(J);
    Scalar vh[9];
    hat3(v, vh);
    Scalar vh2[9];
    mat3_mul(vh, vh, vh2);
    if (theta2 < 1e-8)
    {
        for (int i = 0; i < 9; ++i) J[i] += 0.5 * vh[i] + (1. / 6.) * vh2[i];
        return;
    }
    const Scalar theta = std::sqrt(theta2);
    const Scalar a = (1. - std::cos(theta)) / theta2;
    const Scalar b = (theta - std::sin(theta)) / (theta2 * theta);
    for (int i = 0; i < 9; ++i) J[i] += a * vh[i] + b * vh2[i];
}

// SO(3) right Jacobian J_r(v) = J_l(-v).
inline void so3_right_jacobian(const Scalar v[3], Scalar J[9])
{
    const Scalar nv[3] = {-v[0], -v[1], -v[2]};
    so3_left_jacobian(nv, J);
}

} // namespace examples
} // namespace fatrop

#endif // __fatrop_examples_so3_helpers_hpp__
