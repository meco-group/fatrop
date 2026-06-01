//
// Dense Lie-group example: robot inverse kinematics for a 6-DOF serial arm.
//
// The robot is modelled (with no external dependencies) as a chain of 1-DOF
// revolute joints connected by rigid links. Each joint is described by
//   - a unit rotation axis,                 expressed in its parent link frame,
//   - a location (origin) of the joint,     expressed in its parent link frame,
// and the tool is a fixed offset in the last link frame. This is enough to
// write a custom forward-kinematics expression and its geometric Jacobian by
// hand. The concrete numbers below are the kinematic parameters of a 20 kg /
// 1750 mm-reach collaborative arm (link lengths d1 = 0.236, a2 = 0.862,
// a3 = 0.729, d4 = 0.201, d5 = 0.159, d6 = 0.154 m), with the joint axes and
// link offsets pushed into the base frame at the home (all-zero) configuration.
//
// Inverse kinematics: given a desired end-effector pose (R_des, p_des) we look
// for joint angles q that minimise the squared task-space error
//
//     L(q) = 0.5 || e_p(q) ||^2 + 0.5 || e_R(q) ||^2,
//     e_p(q) = p(q) - p_des                          (position error, world)
//     e_R(q) = Log( R_des^T R(q) )                   (orientation error, so(3))
//
// subject to the joint limits  q_lo <= q <= q_hi. The decision variables (the
// joint angles) are Euclidean, so no manifold retraction is needed; the Lie
// group only shows up in the orientation residual.
//
// The subtle part is d e_R / d q. With the BODY angular velocity omega_b
// (R_dot = R [omega_b]_x) the derivative of the logarithm uses the inverse of
// the *right* Jacobian of SO(3):
//
//     d/dt Log(R_des^T R) = J_r^{-1}( e_R ) * omega_b.
//
// (The target rotation R_des drops out because the body angular velocity of
// R_des^T R equals that of R.) Hence we deliberately build the BODY angular
// Jacobian  omega_b,i = R(q)^T * omega_i^world  so that the *right* inverse
// Jacobian is the correct one. Using the world (left) angular velocity instead
// would require J_l^{-1}.
//
// The cost is a sum of squares, so we use a Gauss-Newton Hessian
// H = J^T J and gradient g = J^T r with the stacked residual r = [e_p; e_R].
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

// ---------------------------------------------------------------------------
// Configurable robot model: a chain of 1-DOF revolute joints.
// ---------------------------------------------------------------------------
struct RevoluteJoint
{
    Scalar axis[3];   // unit rotation axis, in the parent link frame
    Scalar origin[3]; // joint location, in the parent link frame
    Scalar lower;     // joint angle lower limit [rad]
    Scalar upper;     // joint angle upper limit [rad]
};

struct RobotModel
{
    std::vector<RevoluteJoint> joints;
    Scalar tool[3]; // end-effector offset in the last link frame

    Index nq() const { return static_cast<Index>(joints.size()); }
};

// 6-DOF arm with the kinematic parameters of a 20 kg-payload collaborative arm.
// Axes / origins are given in the base frame at the home (q = 0) configuration;
// because the inter-link frames carry no fixed rotation at home, an axis
// expressed in the base equals the same axis expressed in its parent link.
RobotModel make_arm_model()
{
    constexpr Scalar pi = M_PI;
    const Scalar lim = 2.0 * pi; // +/- 360 deg working range on every joint
    RobotModel m;
    m.joints = {
        // axis            origin (location in parent link)      limits
        {{0., 0., 1.}, {0., 0., 0.}, -lim, lim},      // joint 1 (base yaw)
        {{0., -1., 0.}, {0., 0., 0.236}, -lim, lim},  // joint 2 (shoulder)
        {{0., -1., 0.}, {-0.862, 0., 0.}, -lim, lim}, // joint 3 (elbow)
        {{0., -1., 0.}, {-0.729, 0., 0.}, -lim, lim}, // joint 4 (wrist 1)
        {{0., 0., -1.}, {0., -0.201, 0.}, -lim, lim}, // joint 5 (wrist 2)
        {{0., -1., 0.}, {0., 0., -0.159}, -lim, lim}, // joint 6 (wrist 3)
    };
    m.tool[0] = 0.;
    m.tool[1] = -0.154;
    m.tool[2] = 0.;
    return m;
}

// ---------------------------------------------------------------------------
// Forward kinematics + geometric Jacobian (custom, no dependencies).
//
// Walks the chain accumulating the base->link transform (R, p). For every
// joint it records the world-frame rotation axis omega_i and a point o_i on
// that axis, which are exactly what the geometric Jacobian columns need:
//     J_v,i = omega_i x (p_ee - o_i)   (linear velocity of the tool origin)
//     J_w,i = omega_i                  (world angular velocity)
// ---------------------------------------------------------------------------
struct FkResult
{
    Scalar R_ee[9];                       // tool orientation
    Scalar p_ee[3];                       // tool position
    std::vector<std::array<Scalar, 3>> w; // world joint axes
    std::vector<std::array<Scalar, 3>> o; // world point on each joint axis
};

FkResult forward_kinematics(const RobotModel &m, const Scalar *q)
{
    FkResult fk;
    fk.w.resize(m.nq());
    fk.o.resize(m.nq());

    Scalar R[9];
    mat3_identity(R);
    Scalar p[3] = {0., 0., 0.};

    for (Index i = 0; i < m.nq(); ++i)
    {
        const RevoluteJoint &j = m.joints[i];
        // Move to the joint origin in the current link frame: p += R * origin.
        Scalar Ro[3];
        mat3_vec(R, j.origin, Ro);
        p[0] += Ro[0];
        p[1] += Ro[1];
        p[2] += Ro[2];
        // World rotation axis and a point on it, evaluated *before* this joint
        // rotates (the rotation itself does not move its own axis).
        Scalar w_world[3];
        mat3_vec(R, j.axis, w_world);
        for (int k = 0; k < 3; ++k)
        {
            fk.w[i][k] = w_world[k];
            fk.o[i][k] = p[k];
        }
        // Apply the joint rotation: R <- R * Rot(axis, q_i).
        Scalar Rj[9];
        axis_angle_to_mat(j.axis, q[i], Rj);
        Scalar Rnew[9];
        mat3_mul(R, Rj, Rnew);
        for (int k = 0; k < 9; ++k) R[k] = Rnew[k];
    }
    // Tool offset in the last link frame.
    Scalar Rt[3];
    mat3_vec(R, m.tool, Rt);
    for (int k = 0; k < 9; ++k) fk.R_ee[k] = R[k];
    fk.p_ee[0] = p[0] + Rt[0];
    fk.p_ee[1] = p[1] + Rt[1];
    fk.p_ee[2] = p[2] + Rt[2];
    return fk;
}

// ---------------------------------------------------------------------------
// Dense NLP wrapper. Decision variables are the joint angles (Euclidean).
// ---------------------------------------------------------------------------
class RobotIkDense : public DenseAbstract
{
public:
    RobotIkDense(RobotModel model, const Scalar R_des[9], const Scalar p_des[3],
                 std::vector<Scalar> q_init)
        : model_(std::move(model)), q_init_(std::move(q_init))
    {
        for (int k = 0; k < 9; ++k) R_des_[k] = R_des[k];
        for (int k = 0; k < 3; ++k) p_des_[k] = p_des[k];
    }

    Index get_nx() const override { return model_.nq(); }
    Index get_ng() const override { return 0; }
    Index get_ng_ineq() const override { return model_.nq(); } // joint limits

    // Stacked residual r = [e_p (3); e_R (3)] and its 6 x nq Jacobian J.
    void residual_and_jacobian(const Scalar *q, Scalar r[6],
                               std::vector<std::array<Scalar, 6>> &J) const
    {
        const Index nq = model_.nq();
        const FkResult fk = forward_kinematics(model_, q);

        // Position residual e_p = p(q) - p_des (world frame).
        for (int k = 0; k < 3; ++k) r[k] = fk.p_ee[k] - p_des_[k];

        // Orientation residual e_R = Log(R_des^T R(q)).
        Scalar Rdt[9];
        mat3_transpose(R_des_, Rdt);
        Scalar R_err[9];
        mat3_mul(Rdt, fk.R_ee, R_err);
        Scalar q_err[4];
        rotmat_to_quat(R_err, q_err);
        Scalar e_R[3];
        quat_log(q_err, e_R);
        r[3] = e_R[0];
        r[4] = e_R[1];
        r[5] = e_R[2];

        // Inverse RIGHT Jacobian of SO(3) at e_R: maps the body angular
        // velocity to the time derivative of the logarithm.
        Scalar Jri[9];
        so3_right_jacobian_inv(e_R, Jri);
        Scalar Rt[9];
        mat3_transpose(fk.R_ee, Rt); // R(q)^T

        J.assign(nq, {0., 0., 0., 0., 0., 0.});
        for (Index i = 0; i < nq; ++i)
        {
            const Scalar wi[3] = {fk.w[i][0], fk.w[i][1], fk.w[i][2]};
            // Linear part: J_v,i = omega_i x (p_ee - o_i).
            const Scalar d[3] = {fk.p_ee[0] - fk.o[i][0], fk.p_ee[1] - fk.o[i][1],
                                 fk.p_ee[2] - fk.o[i][2]};
            Scalar jv[3];
            cross3(wi, d, jv);
            J[i][0] = jv[0];
            J[i][1] = jv[1];
            J[i][2] = jv[2];
            // Angular part: body axis omega_b = R^T omega_i, then J_r^{-1} * omega_b.
            Scalar wb[3];
            mat3_vec(Rt, wi, wb);
            Scalar jr[3];
            mat3_vec(Jri, wb, jr);
            J[i][3] = jr[0];
            J[i][4] = jr[1];
            J[i][5] = jr[2];
        }
    }

    // Gauss-Newton Hessian H = J^T J and gradient g = J^T r, packed into the
    // (nq + 1) x nq Hht layout (top nq rows = H, last row = g).
    Index eval_Hh(const Scalar *objective_scale, const Scalar *x, const Scalar * /*lam*/,
                  MAT *res) override
    {
        const Index nq = model_.nq();
        blasfeo_gese_wrap(nq + 1, nq, 0.0, res, 0, 0);
        Scalar r[6];
        std::vector<std::array<Scalar, 6>> J;
        residual_and_jacobian(x, r, J);
        const Scalar s = objective_scale[0];
        for (Index a = 0; a < nq; ++a)
        {
            Scalar g = 0.;
            for (int k = 0; k < 6; ++k) g += J[a][k] * r[k];
            blasfeo_matel_wrap(res, nq, a) = s * g;
            for (Index b = 0; b < nq; ++b)
            {
                Scalar h = 0.;
                for (int k = 0; k < 6; ++k) h += J[a][k] * J[b][k];
                blasfeo_matel_wrap(res, a, b) = s * h;
            }
        }
        return 0;
    }

    Index eval_grad(const Scalar *objective_scale, const Scalar *x, Scalar *res) override
    {
        const Index nq = model_.nq();
        Scalar r[6];
        std::vector<std::array<Scalar, 6>> J;
        residual_and_jacobian(x, r, J);
        const Scalar s = objective_scale[0];
        for (Index a = 0; a < nq; ++a)
        {
            Scalar g = 0.;
            for (int k = 0; k < 6; ++k) g += J[a][k] * r[k];
            res[a] = s * g;
        }
        return 0;
    }

    Index eval_f(const Scalar *objective_scale, const Scalar *x, Scalar *res) override
    {
        Scalar r[6];
        std::vector<std::array<Scalar, 6>> J;
        residual_and_jacobian(x, r, J);
        Scalar sum_sq = 0.;
        for (int k = 0; k < 6; ++k) sum_sq += r[k] * r[k];
        *res = 0.5 * objective_scale[0] * sum_sq;
        return 0;
    }

    // No equality constraints.
    Index eval_Ggt(const Scalar *, MAT *res) override
    {
        blasfeo_gese_wrap(res->m, res->n, 0.0, res, 0, 0);
        return 0;
    }
    Index eval_g(const Scalar *, Scalar *) override { return 0; }

    // Joint limits as box inequalities: g_ineq(q) = q, identity Jacobian.
    Index eval_Ggt_ineq(const Scalar *x, MAT *res) override
    {
        const Index nq = model_.nq();
        blasfeo_gese_wrap(nq + 1, nq, 0.0, res, 0, 0);
        for (Index i = 0; i < nq; ++i)
        {
            blasfeo_matel_wrap(res, i, i) = 1.0;   // d g_i / d q_i
            blasfeo_matel_wrap(res, nq, i) = x[i]; // residual row g_i = q_i
        }
        return 0;
    }
    Index eval_gineq(const Scalar *x, Scalar *res) override
    {
        for (Index i = 0; i < model_.nq(); ++i) res[i] = x[i];
        return 0;
    }
    Index get_bounds(Scalar *lower, Scalar *upper) const override
    {
        for (Index i = 0; i < model_.nq(); ++i)
        {
            lower[i] = model_.joints[i].lower;
            upper[i] = model_.joints[i].upper;
        }
        return 0;
    }

    Index get_initial(Scalar *x) const override
    {
        for (Index i = 0; i < model_.nq(); ++i) x[i] = q_init_[i];
        return 0;
    }

private:
    RobotModel model_;
    Scalar R_des_[9];
    Scalar p_des_[3];
    std::vector<Scalar> q_init_;
};

int main()
{
    const RobotModel model = make_arm_model();
    const Index nq = model.nq();

    // Pick a feasible target by running FK on a known joint configuration.
    const std::vector<Scalar> q_true = {0.4, -0.9, 1.2, -0.6, 0.8, 1.1};
    const FkResult target = forward_kinematics(model, q_true.data());

    // Start the solver from a perturbed initial guess (within the joint limits).
    std::mt19937 rng(1);
    std::uniform_real_distribution<Scalar> perturb(-0.4, 0.4);
    std::vector<Scalar> q_init(nq);
    for (Index i = 0; i < nq; ++i) q_init[i] = q_true[i] + perturb(rng);

    auto problem =
        std::make_shared<RobotIkDense>(model, target.R_ee, target.p_ee, q_init);

    OptionRegistry options;
    IpAlgBuilder<DenseType> builder(std::make_shared<NlpDense>(problem));
    auto ipalg = builder.with_options_registry(&options).build();

    Timer timer;
    timer.start();
    IpSolverReturnFlag ret = ipalg->optimize();
    const auto elapsed = timer.stop();
    auto data = builder.get_ipdata();

    std::cout << "Elapsed time: " << elapsed << "\n";
    std::cout << data->timing_statistics() << "\n";

    const VecRealView &q_final = data->current_iterate().primal_x();
    std::vector<Scalar> q_sol(nq);
    for (Index i = 0; i < nq; ++i) q_sol[i] = q_final(i);

    // Verify: pose error of the solution, and joint-limit feasibility.
    const FkResult fk = forward_kinematics(model, q_sol.data());
    const Scalar pos_err =
        std::sqrt((fk.p_ee[0] - target.p_ee[0]) * (fk.p_ee[0] - target.p_ee[0]) +
                  (fk.p_ee[1] - target.p_ee[1]) * (fk.p_ee[1] - target.p_ee[1]) +
                  (fk.p_ee[2] - target.p_ee[2]) * (fk.p_ee[2] - target.p_ee[2]));
    Scalar Rdt[9];
    mat3_transpose(target.R_ee, Rdt);
    Scalar R_err[9];
    mat3_mul(Rdt, fk.R_ee, R_err);
    Scalar q_err[4];
    rotmat_to_quat(R_err, q_err);
    Scalar v_err[3];
    quat_log(q_err, v_err);
    const Scalar rot_err =
        std::sqrt(v_err[0] * v_err[0] + v_err[1] * v_err[1] + v_err[2] * v_err[2]);

    bool within_limits = true;
    for (Index i = 0; i < nq; ++i)
        within_limits =
            within_limits && q_sol[i] >= model.joints[i].lower - 1e-9 &&
            q_sol[i] <= model.joints[i].upper + 1e-9;

    std::cout << "Return flag: " << int(ret) << "\n";
    std::cout << "Success:     " << (ret == IpSolverReturnFlag::Success) << "\n";
    std::cout << "q_solution = [";
    for (Index i = 0; i < nq; ++i)
        std::cout << q_sol[i] << (i + 1 < nq ? ", " : "]\n");
    std::cout << "target position    = [" << target.p_ee[0] << ", " << target.p_ee[1] << ", "
              << target.p_ee[2] << "]\n";
    std::cout << "achieved position  = [" << fk.p_ee[0] << ", " << fk.p_ee[1] << ", "
              << fk.p_ee[2] << "]\n";
    std::cout << "position error     = " << pos_err << " m\n";
    std::cout << "orientation error  = " << rot_err << " rad (~ " << rot_err * 180.0 / M_PI
              << " deg)\n";
    std::cout << "joint limits ok    = " << within_limits << "\n";

    return (ret == IpSolverReturnFlag::Success && pos_err < 1e-6 && rot_err < 1e-6 &&
            within_limits)
               ? 0
               : 1;
}
