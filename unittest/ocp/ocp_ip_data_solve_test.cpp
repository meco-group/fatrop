#include "../random_matrix.hpp"
#include "fatrop/ip_algorithm/ip_data.hpp"
#include "fatrop/linear_algebra/linear_algebra.hpp"
#include "fatrop/nlp/dims.hpp"
#include "fatrop/ocp/aug_system_solver.hpp"
#include "fatrop/ocp/hessian.hpp"
#include "fatrop/ocp/jacobian.hpp"
#include "fatrop/ocp/nlp_ocp.hpp"
#include "fatrop/ocp/ocp_abstract.hpp"
#include "fatrop/ocp/pd_solver_orig.hpp"
#include "fatrop/ocp/pd_system_orig.hpp"
#include <gtest/gtest.h>
#include <memory>

using namespace fatrop;

// example problem 2D point mass
// states: [x, y, vx, vy]
// inputs: [fx, fy]
// dynamics: [xk+1 = xk + dt*vxk, yk+1 = yk + dt*vyk, vxk+1 = vxk + dt*fx/m, vyk+1 = vyk + dt*fy/m]
// cost: fx^2 + fy^2
// constraints:
//  at k = 0: x = 0, y = 0, vx = 0, vy = 0
//  at k = K: x = 1, y = 1, vx = 0, vy = 0
class OcpTest : public OcpAbstract
{
public:
    virtual Index get_nx(const Index k) const { return 4; }

    virtual Index get_nu(const Index k) const
    {
        if (k == 0)
        {
            return 2;
        }
        else if (k == K_ - 1)
        {
            return 0;
        }
        else
        {
            return 2;
        }
    }
    virtual Index get_ng(const Index k) const
    {
        if (k == 0)
        {
            return 4;
        }
        else if (k == K_ - 1)
        {
            return 4;
        }
        else
        {
            return 0;
        }
    };

    virtual Index get_ng_ineq(const Index k) const { return k == K_ - 1 ? 0 : 1; };
    virtual Index get_horizon_length() const { return K_; };
    virtual Index eval_BAbt(const Scalar *states_kp1, const Scalar *inputs_k,
                            const Scalar *states_k, MAT *res, const Index k)
    {
        // set zero
        blasfeo_gese_wrap(res->m, res->n, 0.0, res, 0, 0);
        // Matrix B
        // [  0,    0  ]
        // [  0,    0  ]
        // [ dt/m,  0  ]
        // [  0,   dt/m ]

        // Matrix A
        // [ 1,  0,  dt,  0  ]
        // [ 0,  1,   0,  dt ]
        // [ 0,  0,   1,  0  ]
        // [ 0,  0,   0,  1  ]

        blasfeo_matel_wrap(res, 0, 2) = dt_ / m_;
        blasfeo_matel_wrap(res, 1, 3) = dt_ / m_;

        blasfeo_diare_wrap(4, 1.0, res, 2, 0);
        blasfeo_matel_wrap(res, 4, 0) = dt_;
        blasfeo_matel_wrap(res, 5, 1) = dt_;
        return 0;
    }
    virtual Index eval_RSQrqt(const Scalar *objective_scale, const Scalar *inputs_k,
                              const Scalar *states_k, const Scalar *lam_dyn_k,
                              const Scalar *lam_eq_k, const Scalar *lam_eq_ineq_k, MAT *res,
                              const Index k)
    {
        // set zero
        blasfeo_gese_wrap(res->m, res->n, 0.0, res, 0, 0);
        // Matrix R
        // [ 2,  0 ]
        // [ 0,  2 ]
        if (k < K_ - 1)
        {
            blasfeo_diare_wrap(2, 2.0, res, 0, 0);
        }
        return 0;
    };
    virtual Index eval_Ggt(const Scalar *inputs_k, const Scalar *states_k, MAT *res, const Index k)
    {
        // set zero
        blasfeo_gese_wrap(res->m, res->n, 0.0, res, 0, 0);
        if (k == 0)
            blasfeo_diare_wrap(4, 1.0, res, 2, 0);
        if (k == K_ - 1)
            blasfeo_diare_wrap(4, 1.0, res, 0, 0);
        return 0;
    }
    virtual Index eval_Ggt_ineq(const Scalar *inputs_k, const Scalar *states_k, MAT *res,
                                const Index k)
    {
        if (k == K_ - 1)
            return 0;
        // set zero
        blasfeo_gese_wrap(res->m, res->n, 0.0, res, 0, 0);
        blasfeo_matel_wrap(res, 0, 0) = 1.0;
        return 0;
    };
    virtual Index eval_b(const Scalar *states_kp1, const Scalar *inputs_k, const Scalar *states_k,
                         Scalar *res, const Index k)
    {
        res[0] = -states_kp1[0] + states_k[0] + dt_ * states_k[2] + 0.01 * k;
        res[1] = -states_kp1[1] + states_k[1] + dt_ * states_k[3] + 0.02 * k;
        res[2] = -states_kp1[2] + states_k[2] + dt_ * inputs_k[0] / m_ + 0.03 * k;
        res[3] = -states_kp1[3] + states_k[3] + dt_ * inputs_k[1] / m_ + 0.04 * k;
        return 0;
    }

    virtual Index eval_g(const Scalar *inputs_k, const Scalar *states_k, Scalar *res, const Index k)
    {
        if (k == 0)
        {
            res[0] = states_k[0] + 0.04;
            res[1] = states_k[1] + 0.02;
            res[2] = states_k[2] + 0.01;
            res[3] = states_k[3];
        }
        else if (k == K_ - 1)
        {
            res[0] = states_k[0] - 1.0;
            res[1] = states_k[1] - 2.0;
            res[2] = states_k[2] - 3.0;
            res[3] = states_k[3] - 4.0;
        }
        return 0;
    };
    virtual Index eval_gineq(const Scalar *inputs_k, const Scalar *states_k, Scalar *res,
                             const Index k)
    {
        if (k == K_ - 1)
            return 0;
        res[0] = inputs_k[0];
        return 0;
    };
    virtual Index eval_rq(const Scalar *objective_scale, const Scalar *inputs_k,
                          const Scalar *states_k, Scalar *res, const Index k)
    {
        if (k == K_ - 1)
        {
            res[0] = 0.;
            res[1] = 0.;
            res[2] = 0.;
            res[3] = 0.;
        }
        else
        {
            res[0] = 2 * objective_scale[0] * inputs_k[0] + 0.05 * k;
            res[1] = 2 * objective_scale[0] * inputs_k[1] + 0.06 * k;
            res[2] = 0. + 0.07 * k;
            res[3] = 0. + 0.08 * k;
            res[4] = 0. + 0.09 * k;
            res[5] = 0. + 0.10 * k;
            res[6] = 0. + 0.11 * k;
        }
        return 0;
    }
    virtual Index eval_L(const Scalar *objective_scale, const Scalar *inputs_k,
                         const Scalar *states_k, Scalar *res, const Index k)
    {
        if (k == K_ - 1)
        {
            *res = 0.;
        }
        else
        {
            *res = objective_scale[0] * (inputs_k[0] * inputs_k[0] + inputs_k[1] * inputs_k[1]);
        }
        return 0;
    }
    virtual Index get_bounds(Scalar *lower, Scalar *upper, const Index k) const
    {
        if (k == K_ - 1)
            return 0;
        lower[0] = -0.5;
        upper[0] = 0.5;
        return 0;
    }

    virtual Index get_initial_xk(Scalar *xk, const Index k) const
    {
        xk[0] = 0.0;
        xk[1] = 0.0;
        xk[2] = 0.0;
        xk[3] = 0.0;
        return 0;
    };
    virtual Index get_initial_uk(Scalar *uk, const Index k) const
    {
        if (k == K_ - 1)
            return 0;
        uk[0] = 0.0;
        uk[1] = 0.0;
        return 0;
    };
    virtual ~OcpTest() = default;

private:
    const Index K_ = 3;
    const Scalar m_ = 1.0;
    const Scalar dt_ = 0.05;
};

class OcpImplExampleTest : public ::testing::Test
{
protected:
    OcpImplExampleTest()
        : ocp(std::make_shared<OcpTest>()), nlp(std::make_shared<NlpOcp>(ocp)),
          info(nlp->problem_dims()), data(nlp), D_x(info.number_of_primal_variables),
          D_eq(info.number_of_g_eq_path), D_i(info.number_of_slack_variables),
          aug_solver(std::make_shared<OcpAugSystemSolver>(info)), solver(info, aug_solver),
          rhs_x(info.number_of_primal_variables), rhs_s(info.number_of_slack_variables),
          rhs_g(info.number_of_eq_constraints), rhs_cl(info.number_of_slack_variables),
          rhs_cu(info.number_of_slack_variables)
    {
        data.current_iterate().dual_bounds_l() = 1.;
        data.current_iterate().dual_bounds_u() = 1.;
    }

    std::shared_ptr<OcpTest> ocp;
    std::shared_ptr<NlpOcp> nlp;
    ProblemInfo<OcpType> info;
    IpData<OcpType> data;
    VecRealAllocated D_x, D_eq, D_i;
    std::shared_ptr<OcpAugSystemSolver> aug_solver;
    PdSolverOrig<OcpType> solver;
    VecRealAllocated rhs_x, rhs_s, rhs_g, rhs_cl, rhs_cu;
};

TEST_F(OcpImplExampleTest, SolveLinearSystem)
{
    LinearSystem<PdSystemType<OcpType>> ls(
        info, data.current_iterate().jacobian(), data.current_iterate().hessian(), D_x, false, D_eq,
        D_i, data.current_iterate().delta_lower(), data.current_iterate().delta_upper(),
        data.current_iterate().dual_bounds_l(), data.current_iterate().dual_bounds_u(), rhs_x,
        rhs_s, rhs_g, rhs_cl, rhs_cu);

    EXPECT_NO_THROW(solver.solve_in_place(ls));
}

TEST_F(OcpImplExampleTest, UpdateIterateAndCheckInfeasibility)
{
    rhs_x.block(rhs_x.m(), 0) = data.current_iterate().dual_infeasibility_x();
    rhs_s.block(rhs_s.m(), 0) = data.current_iterate().dual_infeasibility_s();
    rhs_g.block(rhs_g.m(), 0) = data.current_iterate().constr_viol();
    rhs_cl.block(rhs_cl.m(), 0) =
        data.current_iterate().delta_lower() * data.current_iterate().dual_bounds_l();
    rhs_cu.block(rhs_cu.m(), 0) =
        data.current_iterate().delta_upper() * data.current_iterate().dual_bounds_u();
    LinearSystem<PdSystemType<OcpType>> ls(
        info, data.current_iterate().jacobian(), data.current_iterate().hessian(), D_x, false, D_eq,
        D_i, data.current_iterate().delta_lower(), data.current_iterate().delta_upper(),
        data.current_iterate().dual_bounds_l(), data.current_iterate().dual_bounds_u(), rhs_x,
        rhs_s, rhs_g, rhs_cl, rhs_cu);

    LinsolReturnFlag ret = solver.solve_in_place(ls);
    EXPECT_EQ(ret, LinsolReturnFlag::SUCCESS);

    Scalar alpha = 1.0;
    Scalar alpha_z = 1.0;
    data.trial_iterate().primal_x() = data.current_iterate().primal_x() + alpha * rhs_x;
    data.trial_iterate().primal_s() = data.current_iterate().primal_s() + alpha * rhs_s;
    data.trial_iterate().dual_eq() = data.current_iterate().dual_eq() + alpha * rhs_g;
    data.trial_iterate().dual_bounds_l() =
        data.current_iterate().dual_bounds_l() + alpha_z * rhs_cl;
    data.trial_iterate().dual_bounds_u() =
        data.current_iterate().dual_bounds_u() + alpha_z * rhs_cu;

    EXPECT_LT(norm_inf(data.trial_iterate().dual_infeasibility_x()), 1e-6);
    EXPECT_LT(norm_inf(data.trial_iterate().dual_infeasibility_s()), 1e-6);
    EXPECT_LT(norm_inf(data.trial_iterate().constr_viol()), 1e-6);
    // the complementarity constraints are nonlinear so we test the linearized version
    EXPECT_LT(
        norm_inf(data.current_iterate().delta_lower() * data.current_iterate().dual_bounds_l() +
                 data.current_iterate().delta_lower() * rhs_cl +
                 data.current_iterate().dual_bounds_l() * rhs_s),
        1e-6);
    EXPECT_LT(
        norm_inf(data.current_iterate().delta_upper() * data.current_iterate().dual_bounds_u() +
                 data.current_iterate().delta_upper() * rhs_cu +
                 -1. * data.current_iterate().dual_bounds_u() * rhs_s),
        1e-6);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
