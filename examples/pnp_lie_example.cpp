//
// Lie-group OCP example: pinhole-camera PnP (Perspective-n-Point).
//
// We have a known 3D point cloud {P_i} in the world frame and the corresponding 2D
// pinhole projections {u_i} (in normalized image coordinates) observed by a camera with
// unknown pose. We recover the pose T = (R, t) — which maps world points into the
// camera frame as p_cam = R P + t — by minimizing the reprojection error
//
//     L(T) = 0.5 * sum_i || pi(R P_i + t) - u_i ||^2,    pi(X, Y, Z) = (X/Z, Y/Z).
//
// The state x = [t (3), q (4)] lives in R^3 x S^3 (nx = 7); the search direction lives
// in the tangent space [rho (3), phi (3)] = R^6 (nx_tangent = 6). The retraction is a
// product-of-groups (SO(3) x R^3) retraction:
//     t  <- t + alpha * rho               (translation in camera frame)
//     q  <- q (x) exp_quat(0.5 * alpha * phi)  (rotation, right-multiplication)
//
// Linearizing the projection around the current pose:
//     d p_cam = rho - R [P]_x phi,    d r_i = (d pi / d p_cam) * d p_cam,
// where (d pi / d p_cam) is the 2x3 pinhole-projection Jacobian. The Gauss-Newton
// Hessian and gradient are the usual J^T J and J^T r sums.

#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <random>
#include <vector>
#include <fatrop/fatrop.hpp>

using namespace fatrop;

namespace
{
// Hamilton quaternion product: r = a (x) b, both [w, x, y, z].
inline void quat_mul(const Scalar *a, const Scalar *b, Scalar *r)
{
    const Scalar aw = a[0], ax = a[1], ay = a[2], az = a[3];
    const Scalar bw = b[0], bx = b[1], by = b[2], bz = b[3];
    r[0] = aw * bw - ax * bx - ay * by - az * bz;
    r[1] = aw * bx + ax * bw + ay * bz - az * by;
    r[2] = aw * by - ax * bz + ay * bw + az * bx;
    r[3] = aw * bz + ax * by - ay * bx + az * bw;
}

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

inline void quat_exp_half(const Scalar *v, Scalar *q)
{
    const Scalar half_angle = 0.5 * std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (half_angle < 1e-12)
    {
        q[0] = 1. - 0.5 * half_angle * half_angle;
        q[1] = 0.5 * v[0];
        q[2] = 0.5 * v[1];
        q[3] = 0.5 * v[2];
        quat_normalize(q);
        return;
    }
    const Scalar s = std::sin(half_angle) / (2. * half_angle);
    q[0] = std::cos(half_angle);
    q[1] = s * v[0];
    q[2] = s * v[1];
    q[3] = s * v[2];
}

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

// Rotate a 3-vector p by a unit quaternion q.
inline void quat_rotate(const Scalar *q, const Scalar *p, Scalar *p_rot)
{
    const Scalar w = q[0], x = q[1], y = q[2], z = q[3];
    const Scalar vp = x * p[0] + y * p[1] + z * p[2];
    const Scalar vv = x * x + y * y + z * z;
    const Scalar c = 1.0 - 2.0 * vv;
    const Scalar cx = y * p[2] - z * p[1];
    const Scalar cy = z * p[0] - x * p[2];
    const Scalar cz = x * p[1] - y * p[0];
    p_rot[0] = c * p[0] + 2.0 * vp * x + 2.0 * w * cx;
    p_rot[1] = c * p[1] + 2.0 * vp * y + 2.0 * w * cy;
    p_rot[2] = c * p[2] + 2.0 * vp * z + 2.0 * w * cz;
}

inline void cross3(const Scalar *a, const Scalar *b, Scalar *r)
{
    r[0] = a[1] * b[2] - a[2] * b[1];
    r[1] = a[2] * b[0] - a[0] * b[2];
    r[2] = a[0] * b[1] - a[1] * b[0];
}
} // namespace

class PnpLieOcp : public OcpAbstract
{
public:
    PnpLieOcp(std::vector<std::array<Scalar, 3>> world_points,
              std::vector<std::array<Scalar, 2>> image_points,
              const Scalar (&pose_init)[7])
        : world_(std::move(world_points)), image_(std::move(image_points))
    {
        for (int i = 0; i < 7; ++i)
            pose_init_[i] = pose_init[i];
    }

    // State layout: [tx, ty, tz, qw, qx, qy, qz]; tangent: [rho (3), phi (3)].
    Index get_nx(const Index /*k*/) const override { return 7; }
    Index get_nu(const Index /*k*/) const override { return 0; }
    Index get_ng(const Index /*k*/) const override { return 0; }
    Index get_ng_ineq(const Index /*k*/) const override { return 0; }
    Index get_horizon_length() const override { return 1; }

    Index get_nx_tangent(const Index /*k*/) const override { return 6; }
    Index get_nu_tangent(const Index /*k*/) const override { return 0; }

    void apply_retraction_xk(const Index /*k*/, const Scalar *xk, const Scalar *delta_xk,
                             const Scalar alpha, Scalar *xk_next) override
    {
        // Translation: t <- t + alpha * rho.
        xk_next[0] = xk[0] + alpha * delta_xk[0];
        xk_next[1] = xk[1] + alpha * delta_xk[1];
        xk_next[2] = xk[2] + alpha * delta_xk[2];
        // Rotation: q <- q (x) exp_quat(0.5 * alpha * phi).
        const Scalar phi[3] = {alpha * delta_xk[3], alpha * delta_xk[4], alpha * delta_xk[5]};
        Scalar dq[4];
        quat_exp_half(phi, dq);
        quat_mul(xk + 3, dq, xk_next + 3);
        quat_normalize(xk_next + 3);
    }

    Index eval_BAbt(const Scalar *, const Scalar *, const Scalar *, MAT *, const Index) override
    {
        return 0;
    }

    Index eval_RSQrqt(const Scalar *objective_scale, const Scalar * /*inputs_k*/,
                      const Scalar *states_k, const Scalar * /*lam_dyn_k*/,
                      const Scalar * /*lam_eq_k*/, const Scalar * /*lam_eq_ineq_k*/, MAT *res,
                      const Index /*k*/) override
    {
        // Layout: (nu + nx_tan + 1) x (nu + nx_tan) = 7 x 6; top 6x6 = Hessian, last row = grad.
        blasfeo_gese_wrap(7, 6, 0.0, res, 0, 0);

        Scalar H[6][6] = {};
        Scalar g[6] = {};
        accumulate_grad_hess(states_k, g, H);

        const Scalar s = objective_scale[0];
        for (int r = 0; r < 6; ++r)
            for (int c = 0; c < 6; ++c)
                blasfeo_matel_wrap(res, r, c) = s * H[r][c];
        for (int c = 0; c < 6; ++c)
            blasfeo_matel_wrap(res, 6, c) = s * g[c];
        return 0;
    }

    Index eval_Ggt(const Scalar *, const Scalar *, MAT *res, const Index) override
    {
        blasfeo_gese_wrap(res->m, res->n, 0.0, res, 0, 0);
        return 0;
    }
    Index eval_Ggt_ineq(const Scalar *, const Scalar *, MAT *res, const Index) override
    {
        blasfeo_gese_wrap(res->m, res->n, 0.0, res, 0, 0);
        return 0;
    }
    Index eval_b(const Scalar *, const Scalar *, const Scalar *, Scalar *, const Index) override
    {
        return 0;
    }
    Index eval_g(const Scalar *, const Scalar *, Scalar *, const Index) override { return 0; }
    Index eval_gineq(const Scalar *, const Scalar *, Scalar *, const Index) override { return 0; }

    Index eval_rq(const Scalar *objective_scale, const Scalar * /*inputs_k*/,
                  const Scalar *states_k, Scalar *res, const Index /*k*/) override
    {
        Scalar H_unused[6][6] = {};
        Scalar g[6] = {};
        accumulate_grad_hess(states_k, g, H_unused);
        const Scalar s = objective_scale[0];
        for (int i = 0; i < 6; ++i)
            res[i] = s * g[i];
        return 0;
    }

    Index eval_L(const Scalar *objective_scale, const Scalar * /*inputs_k*/,
                 const Scalar *states_k, Scalar *res, const Index /*k*/) override
    {
        const Scalar *t = states_k;
        const Scalar *q = states_k + 3;
        Scalar sum_sq = 0.;
        for (std::size_t i = 0; i < world_.size(); ++i)
        {
            Scalar Rp[3];
            quat_rotate(q, world_[i].data(), Rp);
            const Scalar pcam[3] = {Rp[0] + t[0], Rp[1] + t[1], Rp[2] + t[2]};
            const Scalar inv_z = 1.0 / pcam[2];
            const Scalar du = pcam[0] * inv_z - image_[i][0];
            const Scalar dv = pcam[1] * inv_z - image_[i][1];
            sum_sq += du * du + dv * dv;
        }
        *res = 0.5 * objective_scale[0] * sum_sq;
        return 0;
    }

    Index get_bounds(Scalar *, Scalar *, const Index) const override { return 0; }

    Index get_initial_xk(Scalar *xk, const Index /*k*/) const override
    {
        for (int i = 0; i < 7; ++i)
            xk[i] = pose_init_[i];
        return 0;
    }
    Index get_initial_uk(Scalar *, const Index) const override { return 0; }

private:
    // Builds the 2x6 per-point Jacobian J_i = d r_i / d v with v = [rho; phi], plus residual.
    // p_cam     = R P + t
    // d p_cam   = rho - R [P]_x phi
    // d pi      = (1/Z) [ 1, 0, -X/Z; 0, 1, -Y/Z ]
    // J_i       = d pi * [ I_3 | -R [P]_x ]
    void residual_and_jacobian(const Scalar *t, const Scalar *q, const Scalar *P,
                               const Scalar *u_obs, Scalar r[2], Scalar J[12]) const
    {
        Scalar Rp[3];
        quat_rotate(q, P, Rp);
        const Scalar pcam[3] = {Rp[0] + t[0], Rp[1] + t[1], Rp[2] + t[2]};
        const Scalar inv_z = 1.0 / pcam[2];
        const Scalar inv_z2 = inv_z * inv_z;
        r[0] = pcam[0] * inv_z - u_obs[0];
        r[1] = pcam[1] * inv_z - u_obs[1];

        const Scalar dpi[2][3] = {{inv_z, 0., -pcam[0] * inv_z2},
                                  {0., inv_z, -pcam[1] * inv_z2}};

        // Columns 0..2: d r / d rho = d pi (rho enters p_cam as +I).
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 3; ++j)
                J[6 * i + j] = dpi[i][j];

        // Columns 3..5: d r / d phi = d pi * (-R [P]_x).
        // Column j of (-R [P]_x) = R * (e_j x P), so we rotate the cross product.
        for (int j = 0; j < 3; ++j)
        {
            Scalar e[3] = {0., 0., 0.};
            e[j] = 1.;
            Scalar e_x_P[3];
            cross3(e, P, e_x_P);
            Scalar col[3];
            quat_rotate(q, e_x_P, col);
            for (int i = 0; i < 2; ++i)
                J[6 * i + 3 + j] = dpi[i][0] * col[0] + dpi[i][1] * col[1] + dpi[i][2] * col[2];
        }
    }

    void accumulate_grad_hess(const Scalar *states_k, Scalar g[6], Scalar H[6][6]) const
    {
        const Scalar *t = states_k;
        const Scalar *q = states_k + 3;
        for (std::size_t i = 0; i < world_.size(); ++i)
        {
            Scalar r[2];
            Scalar J[12]; // 2x6 row-major
            residual_and_jacobian(t, q, world_[i].data(), image_[i].data(), r, J);
            for (int row = 0; row < 6; ++row)
            {
                g[row] += J[0 * 6 + row] * r[0] + J[1 * 6 + row] * r[1];
                for (int col = 0; col < 6; ++col)
                    H[row][col] += J[0 * 6 + row] * J[0 * 6 + col] +
                                   J[1 * 6 + row] * J[1 * 6 + col];
            }
        }
    }

    std::vector<std::array<Scalar, 3>> world_;
    std::vector<std::array<Scalar, 2>> image_;
    Scalar pose_init_[7];
};

int main()
{
    // ---- Build a 3D point cloud in the world frame --------------------------------------
    std::mt19937 rng(7);
    std::uniform_real_distribution<Scalar> uniform(-1.0, 1.0);
    constexpr Scalar pixel_noise = 1e-3; // normalized image coords (~1px @ f=1000)
    std::normal_distribution<Scalar> noise(0.0, pixel_noise);

    constexpr int N = 80;
    std::vector<std::array<Scalar, 3>> world(N);
    for (int i = 0; i < N; ++i)
        world[i] = {uniform(rng), uniform(rng), uniform(rng)};

    // ---- True camera pose: world -> camera transform (R, t) ----------------------------
    // True translation places the camera ~5 units back along world -z.
    // True rotation: 20 degrees around axis (0, 1, 0.3)/||.||
    const Scalar axis_raw[3] = {0., 1., 0.3};
    const Scalar axis_n = std::sqrt(axis_raw[0] * axis_raw[0] + axis_raw[1] * axis_raw[1] +
                                    axis_raw[2] * axis_raw[2]);
    const Scalar angle = 20.0 * M_PI / 180.0;
    const Scalar half = 0.5 * angle;
    const Scalar sh = std::sin(half) / axis_n;
    const Scalar q_true[4] = {std::cos(half), axis_raw[0] * sh, axis_raw[1] * sh,
                              axis_raw[2] * sh};
    const Scalar t_true[3] = {0.1, -0.2, 5.0};

    // ---- Project each world point and add image-plane noise ----------------------------
    std::vector<std::array<Scalar, 2>> image(N);
    for (int i = 0; i < N; ++i)
    {
        Scalar Rp[3];
        quat_rotate(q_true, world[i].data(), Rp);
        const Scalar pcam[3] = {Rp[0] + t_true[0], Rp[1] + t_true[1], Rp[2] + t_true[2]};
        const Scalar inv_z = 1.0 / pcam[2];
        image[i] = {pcam[0] * inv_z + noise(rng), pcam[1] * inv_z + noise(rng)};
    }

    // ---- Solve for the camera pose -----------------------------------------------------
    // Initial guess: identity rotation, translation roughly in the ballpark.
    const Scalar pose_init[7] = {0., 0., 4., 1., 0., 0., 0.};

    OptionRegistry options;
    auto ocp = std::make_shared<PnpLieOcp>(world, image, pose_init);
    IpAlgBuilder<OcpType> builder(std::make_shared<NlpOcp>(ocp));
    auto ipalg = builder.with_options_registry(&options).build();

    Timer timer;
    timer.start();
    IpSolverReturnFlag ret = ipalg->optimize();
    const auto elapsed = timer.stop();
    auto data = builder.get_ipdata();
    std::cout << "Elapsed time: " << elapsed << std::endl;
    std::cout << data->timing_statistics() << std::endl;

    const VecRealView &x_final = data->current_iterate().primal_x();
    const Scalar t_est[3] = {x_final(0), x_final(1), x_final(2)};
    const Scalar q_est[4] = {x_final(3), x_final(4), x_final(5), x_final(6)};

    // Recover camera position in world frame: c = -R^T t.
    const Scalar q_inv[4] = {q_est[0], -q_est[1], -q_est[2], -q_est[3]};
    Scalar cam_pos[3];
    {
        const Scalar neg_t[3] = {-t_est[0], -t_est[1], -t_est[2]};
        quat_rotate(q_inv, neg_t, cam_pos);
    }
    const Scalar q_true_inv[4] = {q_true[0], -q_true[1], -q_true[2], -q_true[3]};
    Scalar cam_pos_true[3];
    {
        const Scalar neg_t[3] = {-t_true[0], -t_true[1], -t_true[2]};
        quat_rotate(q_true_inv, neg_t, cam_pos_true);
    }

    std::cout << "Return flag: " << int(ret) << "\n";
    std::cout << "Success:     " << (ret == IpSolverReturnFlag::Success) << "\n";
    std::cout << "t_true  = [" << t_true[0] << ", " << t_true[1] << ", " << t_true[2] << "]\n";
    std::cout << "t_est   = [" << t_est[0] << ", " << t_est[1] << ", " << t_est[2] << "]\n";
    std::cout << "q_true  = [" << q_true[0] << ", " << q_true[1] << ", " << q_true[2] << ", "
              << q_true[3] << "]\n";
    std::cout << "q_est   = [" << q_est[0] << ", " << q_est[1] << ", " << q_est[2] << ", "
              << q_est[3] << "]\n";
    std::cout << "camera position (world frame): true = [" << cam_pos_true[0] << ", "
              << cam_pos_true[1] << ", " << cam_pos_true[2] << "],  est = [" << cam_pos[0]
              << ", " << cam_pos[1] << ", " << cam_pos[2] << "]\n";

    // Geodesic rotation error.
    Scalar q_err[4];
    quat_mul(q_true_inv, q_est, q_err);
    Scalar v_err[3];
    quat_log(q_err, v_err);
    const Scalar rot_err =
        std::sqrt(v_err[0] * v_err[0] + v_err[1] * v_err[1] + v_err[2] * v_err[2]);
    const Scalar pos_err = std::sqrt(
        (cam_pos[0] - cam_pos_true[0]) * (cam_pos[0] - cam_pos_true[0]) +
        (cam_pos[1] - cam_pos_true[1]) * (cam_pos[1] - cam_pos_true[1]) +
        (cam_pos[2] - cam_pos_true[2]) * (cam_pos[2] - cam_pos_true[2]));
    std::cout << "rotation error (rad) = " << rot_err << " (~ " << rot_err * 180.0 / M_PI
              << " deg)\n";
    std::cout << "camera position err  = " << pos_err << "\n";

    // Reprojection RMSE at the recovered pose.
    Scalar sum_sq = 0.;
    for (int i = 0; i < N; ++i)
    {
        Scalar Rp[3];
        quat_rotate(q_est, world[i].data(), Rp);
        const Scalar pcam[3] = {Rp[0] + t_est[0], Rp[1] + t_est[1], Rp[2] + t_est[2]};
        const Scalar inv_z = 1.0 / pcam[2];
        const Scalar du = pcam[0] * inv_z - image[i][0];
        const Scalar dv = pcam[1] * inv_z - image[i][1];
        sum_sq += du * du + dv * dv;
    }
    const Scalar rmse = std::sqrt(sum_sq / N);
    std::cout << "reprojection RMSE    = " << rmse << " (pixel noise sigma = " << pixel_noise
              << ")\n";

    return (ret == IpSolverReturnFlag::Success && rot_err < 0.02 && pos_err < 0.05) ? 0 : 1;
}
