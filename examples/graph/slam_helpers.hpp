//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//
// Visual-SLAM math helpers built on top of so3_helpers.hpp:
//
//   * std::array typed aliases for SE(3) poses, landmarks and 3x3 matrices.
//   * SE(3) helpers that return / consume those typed aliases
//     (quat_to_rot, rot_t_apply, rot_to_quat, hat).
//   * Pinhole camera projection and its 2x6/2x3 Jacobian (in tangent space).
//   * Body-frame relative-pose residual and its 6x6/6x6 Jacobian.
//
// Everything lives in fatrop::examples to share the namespace with
// so3_helpers.hpp.
//

#ifndef __fatrop_examples_slam_helpers_hpp__
#define __fatrop_examples_slam_helpers_hpp__

#include "fatrop/context/context.hpp"
#include "so3_helpers.hpp"

#include <array>
#include <cmath>

namespace fatrop
{
namespace examples
{

// ============================================================================
// Type aliases
// ============================================================================
using Vec3 = std::array<Scalar, 3>;
using Quat = std::array<Scalar, 4>;
using Mat3 = std::array<std::array<Scalar, 3>, 3>;

constexpr Index POSE_PRIM_DIM = 7; // (px, py, pz, qw, qx, qy, qz)
constexpr Index POSE_TAN_DIM = 6;  // (delta_p (3) ; delta_phi (3))
constexpr Index LMK_DIM = 3;       // (lx, ly, lz)

using Pose = std::array<Scalar, POSE_PRIM_DIM>;
using Landmark = std::array<Scalar, LMK_DIM>;

// ============================================================================
// SLAM measurement types
// ============================================================================
// Generic — used by any consumer of odometry / pinhole observations
// (synthetic data generators, problem formulations, ...). The semantics:
//
//   OdomMeas:    measured body-frame relative pose between two pose
//                indices a and b. The 6-vector z stacks
//                (delta_p_body (3) ; delta_phi_body (3)) using the same
//                convention as rel_pose_raw().
//   PixelMeas:   measured pinhole projection (u, v) of landmark lmk_j
//                seen from pose pose_k, in the same convention as
//                project().
struct OdomMeas
{
    Index a, b;
    Scalar z[6];
};

struct PixelMeas
{
    Index pose_k;
    Index lmk_j;
    Scalar uv[2];
};

// ============================================================================
// SE(3) helpers — typed wrappers over so3_helpers.hpp
// ============================================================================

// Unit quaternion -> 3x3 rotation matrix (Mat3 layout).
inline Mat3 quat_to_rot(const Scalar *q)
{
    Scalar Rm[9];
    quat_to_rotmat(q, Rm);
    return {{{Rm[0], Rm[1], Rm[2]}, {Rm[3], Rm[4], Rm[5]}, {Rm[6], Rm[7], Rm[8]}}};
}

// Apply R^T (i.e. rotate from world frame into body frame).
inline Vec3 rot_t_apply(const Mat3 &R, const Vec3 &x)
{
    return {R[0][0] * x[0] + R[1][0] * x[1] + R[2][0] * x[2],
            R[0][1] * x[0] + R[1][1] * x[1] + R[2][1] * x[2],
            R[0][2] * x[0] + R[1][2] * x[1] + R[2][2] * x[2]};
}

// q_out = q_a^{-1} * q_b. Convenient for body-frame relative rotations.
inline void quat_inv_mul(const Scalar *a, const Scalar *b, Scalar *out)
{
    const Scalar a_inv[4] = {a[0], -a[1], -a[2], -a[3]};
    quat_mul(a_inv, b, out);
}

// Skew-symmetric "hat" map [v]_x written into Mat3 form.
inline void hat(const Scalar *v, Mat3 &S)
{
    S[0] = {0.0, -v[2], v[1]};
    S[1] = {v[2], 0.0, -v[0]};
    S[2] = {-v[1], v[0], 0.0};
}

// 3x3 rotation matrix -> unit quaternion (Shepperd's method).
inline Quat rot_to_quat(const Mat3 &R)
{
    const Scalar tr = R[0][0] + R[1][1] + R[2][2];
    Quat q;
    if (tr > 0.0)
    {
        const Scalar s = std::sqrt(tr + 1.0) * 2.0;
        q[0] = 0.25 * s;
        q[1] = (R[2][1] - R[1][2]) / s;
        q[2] = (R[0][2] - R[2][0]) / s;
        q[3] = (R[1][0] - R[0][1]) / s;
    }
    else if (R[0][0] > R[1][1] && R[0][0] > R[2][2])
    {
        const Scalar s = std::sqrt(1.0 + R[0][0] - R[1][1] - R[2][2]) * 2.0;
        q[0] = (R[2][1] - R[1][2]) / s;
        q[1] = 0.25 * s;
        q[2] = (R[0][1] + R[1][0]) / s;
        q[3] = (R[0][2] + R[2][0]) / s;
    }
    else if (R[1][1] > R[2][2])
    {
        const Scalar s = std::sqrt(1.0 + R[1][1] - R[0][0] - R[2][2]) * 2.0;
        q[0] = (R[0][2] - R[2][0]) / s;
        q[1] = (R[0][1] + R[1][0]) / s;
        q[2] = 0.25 * s;
        q[3] = (R[1][2] + R[2][1]) / s;
    }
    else
    {
        const Scalar s = std::sqrt(1.0 + R[2][2] - R[0][0] - R[1][1]) * 2.0;
        q[0] = (R[1][0] - R[0][1]) / s;
        q[1] = (R[0][2] + R[2][0]) / s;
        q[2] = (R[1][2] + R[2][1]) / s;
        q[3] = 0.25 * s;
    }
    return q;
}

// Inverse right Jacobian of SO(3) at phi (3x3, row-major into out):
//
//   J_r^{-1}(phi) = I + (1/2)[phi]_x + c(theta) [phi]_x^2 ,
//   c(theta)      = 1/theta^2 - cos(theta/2) / (2 theta sin(theta/2)) ,
//
// with a Taylor branch below |phi|^2 < 1e-8 to avoid the 0/0 form. The
// inverse *left* Jacobian is simply the transpose:
//   J_l^{-1}(phi) = J_r^{-1}(phi)^T .
inline void right_jac_inv_so3(const Scalar phi[3], Scalar out[3][3])
{
    const Scalar theta2 = phi[0] * phi[0] + phi[1] * phi[1] + phi[2] * phi[2];
    Scalar c;
    if (theta2 < 1e-8)
    {
        // Taylor: c(theta) = 1/12 + theta^2/720 + O(theta^4)
        c = (1.0 / 12.0) + theta2 * (1.0 / 720.0);
    }
    else
    {
        const Scalar theta = std::sqrt(theta2);
        const Scalar half = 0.5 * theta;
        c = 1.0 / theta2 - std::cos(half) / (2.0 * theta * std::sin(half));
    }
    for (Index i = 0; i < 3; ++i)
        for (Index j = 0; j < 3; ++j)
            out[i][j] = (i == j) ? 1.0 : 0.0;
    // + (1/2) [phi]_x
    out[0][1] += -0.5 * phi[2];
    out[0][2] += 0.5 * phi[1];
    out[1][0] += 0.5 * phi[2];
    out[1][2] += -0.5 * phi[0];
    out[2][0] += -0.5 * phi[1];
    out[2][1] += 0.5 * phi[0];
    // + c * [phi]_x^2 == c * (phi phi^T - theta^2 I)
    for (Index i = 0; i < 3; ++i)
        for (Index j = 0; j < 3; ++j)
            out[i][j] += c * phi[i] * phi[j];
    for (Index i = 0; i < 3; ++i)
        out[i][i] -= c * theta2;
}

// ============================================================================
// Body-frame relative-pose residual (odometry edges)
// ============================================================================
// r = ( R(q_i)^T (p_j - p_i)               ; "delta_p_body"
//       log( q_i^{-1} q_j )                 ) "delta_phi_body"
// xi, xj point at the 7-Scalar primal (px, py, pz, qw, qx, qy, qz).

inline void rel_pose_raw(const Scalar *xi, const Scalar *xj, Scalar *out)
{
    const Mat3 Ri = quat_to_rot(xi + 3);
    const Vec3 dp{xj[0] - xi[0], xj[1] - xi[1], xj[2] - xi[2]};
    const Vec3 dp_body = rot_t_apply(Ri, dp);
    out[0] = dp_body[0];
    out[1] = dp_body[1];
    out[2] = dp_body[2];
    Scalar q_rel[4];
    quat_inv_mul(xi + 3, xj + 3, q_rel);
    quat_log(q_rel, out + 3);
}

// Jacobian of rel_pose_raw w.r.t. the 6D Euclidean / right-trivialised
// tangent perturbations of xi (Ji) and xj (Jj). 
inline void rel_pose_jacobian(const Scalar *xi, const Scalar *xj, Scalar Ji[6][6],
                              Scalar Jj[6][6])
{
    for (Index i = 0; i < 6; ++i)
        for (Index j = 0; j < 6; ++j)
        {
            Ji[i][j] = 0.0;
            Jj[i][j] = 0.0;
        }
    const Mat3 Ri = quat_to_rot(xi + 3);
    const Vec3 dp{xj[0] - xi[0], xj[1] - xi[1], xj[2] - xi[2]};
    const Vec3 dp_body = rot_t_apply(Ri, dp);
    Mat3 dp_body_hat;
    hat(dp_body.data(), dp_body_hat);

    for (Index r = 0; r < 3; ++r)
    {
        // d r_t / d delta_p_world (left:  -R^T, right: +R^T)
        for (Index c = 0; c < 3; ++c)
        {
            Ji[r][c] = -Ri[c][r];
            Jj[r][c] = Ri[c][r];
        }
        // d r_t / d delta_phi_a = [dp_body]_x
        for (Index c = 0; c < 3; ++c)
            Ji[r][3 + c] = dp_body_hat[r][c];
    }
    // ---- Rotation rows (exact via right/left SO(3) Jacobian inverses) --
    //   r_phi = log( q_i^{-1} q_j )
    //   d r_phi / d delta_phi_b =  J_r^{-1}(r_phi)
    //   d r_phi / d delta_phi_a = -J_l^{-1}(r_phi) = -J_r^{-1}(r_phi)^T
    Scalar q_rel[4];
    quat_inv_mul(xi + 3, xj + 3, q_rel);
    Scalar r_phi[3];
    quat_log(q_rel, r_phi);
    Scalar Jri[3][3];
    right_jac_inv_so3(r_phi, Jri);
    for (Index i = 0; i < 3; ++i)
        for (Index j = 0; j < 3; ++j)
        {
            Jj[3 + i][3 + j] = Jri[i][j];  //  J_r^{-1}(r_phi)
            Ji[3 + i][3 + j] = -Jri[j][i]; // -J_l^{-1}(r_phi) = -J_r^{-1}(r_phi)^T
        }
}

// ============================================================================
// Pinhole camera projection (intrinsic constants below)
// ============================================================================
// Camera coordinates (z_cam, x_cam, y_cam) = R(q_k)^T (l_j - p_k).
// Pixel measurement: u = f * x_cam / z_cam, v = f * y_cam / z_cam.

struct CameraIntrinsics
{
    Scalar focal = 200.0; // pixels / unit
    Scalar z_min = 0.1;   // numerical floor below which a landmark is "behind" the camera
};

// Default intrinsics used by the helpers below. Tweak if the example wants
// a different camera model.
constexpr Scalar FOCAL = 200.0;
constexpr Scalar Z_MIN = 0.1;

// Project a landmark into a pose's camera. Returns false if the point is
// behind the camera (l_cam_z <= Z_MIN).
inline bool project(const Pose &pose, const Landmark &lmk, Scalar pix[2])
{
    const Mat3 R = quat_to_rot(pose.data() + 3);
    const Vec3 dp{lmk[0] - pose[0], lmk[1] - pose[1], lmk[2] - pose[2]};
    const Vec3 lc = rot_t_apply(R, dp);
    if (lc[2] <= Z_MIN)
        return false;
    pix[0] = FOCAL * lc[0] / lc[2];
    pix[1] = FOCAL * lc[1] / lc[2];
    return true;
}

// 2x6 pose Jacobian Jp and 2x3 landmark Jacobian Jl of the pinhole
// projection at (pose, lmk). Returns false if the landmark is not in
// front of the camera (no Jacobian written).
inline bool project_jacobian(const Pose &pose, const Landmark &lmk, Scalar Jp[2][6],
                             Scalar Jl[2][3])
{
    const Mat3 R = quat_to_rot(pose.data() + 3);
    const Vec3 dp{lmk[0] - pose[0], lmk[1] - pose[1], lmk[2] - pose[2]};
    const Vec3 lc = rot_t_apply(R, dp);
    if (lc[2] <= Z_MIN)
        return false;
    const Scalar inv_z = 1.0 / lc[2];
    const Scalar inv_z2 = inv_z * inv_z;
    Mat3 lc_hat;
    hat(lc.data(), lc_hat);
    for (Index a = 0; a < 6; ++a)
    {
        Scalar dlx, dly, dlz;
        if (a < 3)
        {
            // d l_cam / d delta_p_world[a] = -R^T column a
            dlx = -R[a][0];
            dly = -R[a][1];
            dlz = -R[a][2];
        }
        else
        {
            // d l_cam / d delta_phi[a-3] = column (a-3) of [l_cam]_x
            const Index c = a - 3;
            dlx = lc_hat[0][c];
            dly = lc_hat[1][c];
            dlz = lc_hat[2][c];
        }
        Jp[0][a] = FOCAL * (dlx * lc[2] - lc[0] * dlz) * inv_z2;
        Jp[1][a] = FOCAL * (dly * lc[2] - lc[1] * dlz) * inv_z2;
    }
    for (Index a = 0; a < 3; ++a)
    {
        // d l_cam / d landmark[a] = R^T column a
        const Scalar dlx = R[a][0];
        const Scalar dly = R[a][1];
        const Scalar dlz = R[a][2];
        Jl[0][a] = FOCAL * (dlx * lc[2] - lc[0] * dlz) * inv_z2;
        Jl[1][a] = FOCAL * (dly * lc[2] - lc[1] * dlz) * inv_z2;
    }
    return true;
}

} // namespace examples
} // namespace fatrop

#endif // __fatrop_examples_slam_helpers_hpp__
