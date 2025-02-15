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
#include "ocp_test_probem.hpp"
#include <gtest/gtest.h>
#include <memory>

using namespace fatrop;
using namespace fatrop::test;

// example problem 2D point mass
// states: [x, y, vx, vy]
// inputs: [fx, fy]
// dynamics: [xk+1 = xk + dt*vxk, yk+1 = yk + dt*vyk, vxk+1 = vxk + dt*fx/m, vyk+1 = vyk + dt*fy/m]
// cost: fx^2 + fy^2
// constraints:
//  at k = 0: x = 0, y = 0, vx = 0, vy = 0
//  at k = K: x = 1, y = 1, vx = 0, vy = 0
class OcpImplExampleTest : public ::testing::Test
{
protected:
    OcpImplExampleTest()
        : ocp(std::make_shared<OcpTestProblem>()), nlp(std::make_shared<NlpOcp>(ocp)),
          info(nlp->problem_dims()), data(nlp), D_x(info.number_of_primal_variables + info.number_of_slack_variables),
          D_eq(info.number_of_eq_constraints),
          aug_solver(std::make_shared<AugSystemSolver<OcpType>>(info)), solver(info, aug_solver),
          rhs_x(info.number_of_primal_variables), rhs_s(info.number_of_slack_variables),
          rhs_g(info.number_of_eq_constraints), rhs_cl(info.number_of_slack_variables),
          rhs_cu(info.number_of_slack_variables)
    {
        data.current_iterate().set_dual_bounds_l(
            VecRealScalar(data.current_iterate().dual_bounds_l().m(), 1));
        data.current_iterate().set_dual_bounds_u(
            VecRealScalar(data.current_iterate().dual_bounds_l().m(), 1));
    }

    std::shared_ptr<OcpTestProblem> ocp;
    std::shared_ptr<NlpOcp> nlp;
    ProblemInfo<OcpType> info;
    IpData<OcpType> data;
    VecRealAllocated D_x, D_eq;
    std::shared_ptr<AugSystemSolver<OcpType>> aug_solver;
    PdSolverOrig<OcpType> solver;
    VecRealAllocated rhs_x, rhs_s, rhs_g, rhs_cl, rhs_cu;
};

TEST_F(OcpImplExampleTest, SolveLinearSystem)
{
    LinearSystem<PdSystemType<OcpType>> ls(
        info, data.current_iterate().jacobian(), data.current_iterate().hessian(), D_x, true, D_eq,
        data.current_iterate().delta_lower(), data.current_iterate().delta_upper(),
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
        info, data.current_iterate().jacobian(), data.current_iterate().hessian(), D_x, true, D_eq,
        data.current_iterate().delta_lower(), data.current_iterate().delta_upper(),
        data.current_iterate().dual_bounds_l(), data.current_iterate().dual_bounds_u(), rhs_x,
        rhs_s, rhs_g, rhs_cl, rhs_cu);

    LinsolReturnFlag ret = solver.solve_in_place(ls);
    EXPECT_EQ(ret, LinsolReturnFlag::SUCCESS);
    VecRealAllocated x_full(ls.m());
    ls.get_rhs(x_full);
    EXPECT_TRUE(x_full.is_finite());

    Scalar alpha = 1.0;
    Scalar alpha_z = 1.0;
    data.trial_iterate().set_primal_x(data.current_iterate().primal_x() + alpha * rhs_x);
    data.trial_iterate().set_primal_s(data.current_iterate().primal_s() + alpha * rhs_s);
    data.trial_iterate().set_dual_eq(data.current_iterate().dual_eq() + alpha * rhs_g);
    data.trial_iterate().set_dual_bounds_l(data.current_iterate().dual_bounds_l() +
                                           alpha_z * rhs_cl);
    data.trial_iterate().set_dual_bounds_u(data.current_iterate().dual_bounds_u() +
                                           alpha_z * rhs_cu);
    data.store_current_iterate();
    data.accept_trial_iterate();

    EXPECT_LT(norm_inf(data.current_iterate().dual_infeasibility_x()), 1e-6);
    EXPECT_LT(norm_inf(data.current_iterate().dual_infeasibility_s()), 1e-6);
    EXPECT_LT(norm_inf(data.current_iterate().constr_viol()), 1e-6);

    data.restore_current_iterate();
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
