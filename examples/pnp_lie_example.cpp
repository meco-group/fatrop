//
// Dense Lie-group example: pinhole-camera PnP (Perspective-n-Point).
//
// We have a known 3D point cloud {P_i} in the world frame and the corresponding
// 2D pinhole projections {u_i} (in normalized image coordinates) observed by a
// camera with unknown pose. We recover the pose T = (R, t) — which maps world
// points into the camera frame as p_cam = R P + t — by minimizing the
// reprojection error
//
//     L(T) = 0.5 * sum_i || pi(R P_i + t) - u_i ||^2,    pi(X, Y, Z) = (X/Z, Y/Z).
//
// The state x = [t (3), q (4)] lives in R^3 x S^3 (get_nx() == 7); the search
// direction lives in [rho (3), phi (3)] (get_nx_tangent() == 6). The retraction
// is the usual product-of-groups update:
//     t  <- t + alpha * rho
//     q  <- q (x) Exp_q(alpha * phi)
//

#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <random>
#include <vector>
#include <fatrop/fatrop.hpp>

#include "so3_helpers.hpp"

using namespace fatrop;
using namespace fatrop::examples;

class PnpLieDense : public DenseAbstract
{
public:
    PnpLieDense(std::vector<std::array<Scalar, 3>> world_points,
                std::vector<std::array<Scalar, 2>> image_points, const Scalar (&pose_init)[7])
        : world_(std::move(world_points)), image_(std::move(image_points))
    {
        for (int i = 0; i < 7; ++i)
            pose_init_[i] = pose_init[i];
    }

    // State: [tx, ty, tz, qw, qx, qy, qz] (nx = 7).
    // Tangent: [rho (3), phi (3)]                (nx_tangent = 6).
    Index get_nx() const override { return 7; }
    Index get_nx_tangent() const override { return 6; }
    Index get_ng() const override { return 0; }
    Index get_ng_ineq() const override { return 0; }

    void apply_retraction(const Scalar *x, const Scalar *delta, const Scalar alpha,
                          Scalar *x_next) override
    {
        // Translation update.
        x_next[0] = x[0] + alpha * delta[0];
        x_next[1] = x[1] + alpha * delta[1];
        x_next[2] = x[2] + alpha * delta[2];
        // Rotation update: right-multiply by exp_quat(0.5 * alpha * phi).
        const Scalar phi[3] = {alpha * delta[3], alpha * delta[4], alpha * delta[5]};
        Scalar dq[4];
        quat_exp(phi, dq);
        quat_mul(x + 3, dq, x_next + 3);
        quat_normalize(x_next + 3);
    }

    Index eval_Hh(const Scalar *objective_scale, const Scalar *x, const Scalar * /*lam*/,
                  MAT *res) override
    {
        blasfeo_gese_wrap(7, 6, 0.0, res, 0, 0);
        Scalar H[6][6] = {};
        Scalar g[6] = {};
        accumulate_grad_hess(x, g, H);
        const Scalar s = objective_scale[0];
        for (int r = 0; r < 6; ++r)
            for (int c = 0; c < 6; ++c)
                blasfeo_matel_wrap(res, r, c) = s * H[r][c];
        for (int c = 0; c < 6; ++c)
            blasfeo_matel_wrap(res, 6, c) = s * g[c];
        return 0;
    }

    Index eval_Ggt(const Scalar *, MAT *res) override
    {
        blasfeo_gese_wrap(res->m, res->n, 0.0, res, 0, 0);
        return 0;
    }
    Index eval_Ggt_ineq(const Scalar *, MAT *res) override
    {
        blasfeo_gese_wrap(res->m, res->n, 0.0, res, 0, 0);
        return 0;
    }
    Index eval_g(const Scalar *, Scalar *) override { return 0; }
    Index eval_gineq(const Scalar *, Scalar *) override { return 0; }

    Index eval_grad(const Scalar *objective_scale, const Scalar *x, Scalar *res) override
    {
        Scalar H_unused[6][6] = {};
        Scalar g[6] = {};
        accumulate_grad_hess(x, g, H_unused);
        const Scalar s = objective_scale[0];
        for (int i = 0; i < 6; ++i)
            res[i] = s * g[i];
        return 0;
    }

    Index eval_f(const Scalar *objective_scale, const Scalar *x, Scalar *res) override
    {
        const Scalar *t = x;
        const Scalar *q = x + 3;
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

    Index get_bounds(Scalar *, Scalar *) const override { return 0; }
    Index get_initial(Scalar *x) const override
    {
        for (int i = 0; i < 7; ++i)
            x[i] = pose_init_[i];
        return 0;
    }

private:
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
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 3; ++j)
                J[6 * i + j] = dpi[i][j];
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

    void accumulate_grad_hess(const Scalar *x, Scalar g[6], Scalar H[6][6]) const
    {
        const Scalar *t = x;
        const Scalar *q = x + 3;
        for (std::size_t i = 0; i < world_.size(); ++i)
        {
            Scalar r[2];
            Scalar J[12];
            residual_and_jacobian(t, q, world_[i].data(), image_[i].data(), r, J);
            for (int row = 0; row < 6; ++row)
            {
                g[row] += J[0 * 6 + row] * r[0] + J[1 * 6 + row] * r[1];
                for (int col = 0; col < 6; ++col)
                    H[row][col] +=
                        J[0 * 6 + row] * J[0 * 6 + col] + J[1 * 6 + row] * J[1 * 6 + col];
            }
        }
    }

    std::vector<std::array<Scalar, 3>> world_;
    std::vector<std::array<Scalar, 2>> image_;
    Scalar pose_init_[7];
};

int main()
{
    std::mt19937 rng(7);
    std::uniform_real_distribution<Scalar> uniform(-1.0, 1.0);
    constexpr Scalar pixel_noise = 1e-3;
    std::normal_distribution<Scalar> noise(0.0, pixel_noise);

    constexpr int N = 80;
    std::vector<std::array<Scalar, 3>> world(N);
    for (int i = 0; i < N; ++i)
        world[i] = {uniform(rng), uniform(rng), uniform(rng)};

    const Scalar axis_raw[3] = {0., 1., 0.3};
    const Scalar axis_n = std::sqrt(axis_raw[0] * axis_raw[0] + axis_raw[1] * axis_raw[1] +
                                    axis_raw[2] * axis_raw[2]);
    const Scalar angle = 20.0 * M_PI / 180.0;
    const Scalar half = 0.5 * angle;
    const Scalar sh = std::sin(half) / axis_n;
    const Scalar q_true[4] = {std::cos(half), axis_raw[0] * sh, axis_raw[1] * sh,
                              axis_raw[2] * sh};
    const Scalar t_true[3] = {0.1, -0.2, 5.0};

    std::vector<std::array<Scalar, 2>> image(N);
    for (int i = 0; i < N; ++i)
    {
        Scalar Rp[3];
        quat_rotate(q_true, world[i].data(), Rp);
        const Scalar pcam[3] = {Rp[0] + t_true[0], Rp[1] + t_true[1], Rp[2] + t_true[2]};
        const Scalar inv_z = 1.0 / pcam[2];
        image[i] = {pcam[0] * inv_z + noise(rng), pcam[1] * inv_z + noise(rng)};
    }

    const Scalar pose_init[7] = {0., 0., 4., 1., 0., 0., 0.};

    OptionRegistry options;
    auto problem = std::make_shared<PnpLieDense>(world, image, pose_init);
    IpAlgBuilder<DenseType> builder(std::make_shared<NlpDense>(problem));
    auto ipalg = builder.with_options_registry(&options).build();

    Timer timer;
    timer.start();
    IpSolverReturnFlag ret = ipalg->optimize();
    const auto elapsed = timer.stop();
    auto data = builder.get_ipdata();
    std::cout << "Elapsed time: " << elapsed << "\n";
    std::cout << data->timing_statistics() << "\n";

    const VecRealView &x_final = data->current_iterate().primal_x();
    const Scalar t_est[3] = {x_final(0), x_final(1), x_final(2)};
    const Scalar q_est[4] = {x_final(3), x_final(4), x_final(5), x_final(6)};
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

    Scalar q_err[4];
    quat_mul(q_true_inv, q_est, q_err);
    Scalar v_err[3];
    quat_log(q_err, v_err);
    const Scalar rot_err =
        std::sqrt(v_err[0] * v_err[0] + v_err[1] * v_err[1] + v_err[2] * v_err[2]);
    const Scalar pos_err =
        std::sqrt((cam_pos[0] - cam_pos_true[0]) * (cam_pos[0] - cam_pos_true[0]) +
                  (cam_pos[1] - cam_pos_true[1]) * (cam_pos[1] - cam_pos_true[1]) +
                  (cam_pos[2] - cam_pos_true[2]) * (cam_pos[2] - cam_pos_true[2]));
    std::cout << "rotation error (rad) = " << rot_err << " (~ " << rot_err * 180.0 / M_PI
              << " deg)\n";
    std::cout << "camera position err  = " << pos_err << "\n";

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
