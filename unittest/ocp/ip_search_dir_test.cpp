#include "../random_matrix.hpp"
#include "fatrop/ip_algorithm/ip_data.hpp"
#include "fatrop/ip_algorithm/ip_search_dir.hpp"
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
class IpSearchDirTest : public ::testing::Test
{
protected:
    IpSearchDirTest()
        : ocp(std::make_shared<OcpTestProblem>()), nlp(std::make_shared<NlpOcp>(ocp)),
          info(nlp->problem_dims()), data(std::make_shared<IpData<OcpType>>(nlp)),
          D_x(info.number_of_primal_variables), D_eq(info.number_of_g_eq_path),
          D_i(info.number_of_slack_variables),
          aug_solver(std::make_shared<AugSystemSolver<OcpType>>(info)),
          solver(std::make_shared<PdSolverOrig<OcpType>>(info, aug_solver)),
          search_dir(data, solver)
    {
        data->current_iterate().set_mu(1.0);
        data->current_iterate().set_dual_bounds_l(
            VecRealScalar(data->current_iterate().dual_bounds_l().m(), 1));
        data->current_iterate().set_dual_bounds_u(
            VecRealScalar(data->current_iterate().dual_bounds_l().m(), 1));
    }

    std::shared_ptr<OcpTestProblem> ocp;
    std::shared_ptr<NlpOcp> nlp;
    ProblemInfo<OcpType> info;
    std::shared_ptr<IpData<OcpType>> data;
    VecRealAllocated D_x, D_eq, D_i;
    std::shared_ptr<AugSystemSolver<OcpType>> aug_solver;
    std::shared_ptr<PdSolverOrig<OcpType>> solver;
    IpSearchDirImpl<OcpType> search_dir;
};

TEST_F(IpSearchDirTest, SolveLinearSystem) { EXPECT_NO_THROW(search_dir.compute_search_dir()); }

TEST_F(IpSearchDirTest, UpdateIterateAndCheckInfeasibility)
{
    LinsolReturnFlag ret = search_dir.compute_search_dir();
    EXPECT_EQ(ret, LinsolReturnFlag::SUCCESS);

    Scalar alpha = 1.0;
    Scalar alpha_z = 1.0;
    Scalar mu = data->current_iterate().mu();
    data->trial_iterate().set_primal_x(data->current_iterate().primal_x() +
                                       alpha * data->current_iterate().delta_primal_x());
    data->trial_iterate().set_primal_s(data->current_iterate().primal_s() +
                                       alpha * data->current_iterate().delta_primal_s());
    data->trial_iterate().set_dual_eq(data->current_iterate().dual_eq() +
                                      alpha * data->current_iterate().delta_dual_eq());
    data->trial_iterate().set_dual_bounds_l(data->current_iterate().dual_bounds_l() +
                                            alpha_z *
                                                data->current_iterate().delta_dual_bounds_l());
    data->trial_iterate().set_dual_bounds_u(data->current_iterate().dual_bounds_u() +
                                            alpha_z *
                                                data->current_iterate().delta_dual_bounds_u());

    EXPECT_LT(norm_inf(data->trial_iterate().dual_infeasibility_x()), 1e-6);
    EXPECT_LT(norm_inf(data->trial_iterate().dual_infeasibility_s()), 1e-6);
    EXPECT_LT(norm_inf(data->trial_iterate().constr_viol()), 1e-6);
    const Index m = data->current_iterate().primal_s().m();
    // the complementarity constraints are nonlinear so we test the linearized version
    EXPECT_LT(
        norm_inf(
            if_else(data->current_iterate().lower_bounded(),
                    data->current_iterate().delta_lower() * data->current_iterate().dual_bounds_l() - mu,
                    VecRealScalar(m, 0.)) +
            data->current_iterate().delta_lower() * data->current_iterate().delta_dual_bounds_l() +
            data->current_iterate().dual_bounds_l() * data->current_iterate().delta_primal_s()),
        1e-6);
    EXPECT_LT(
        norm_inf(
            if_else(data->current_iterate().upper_bounded(),
                    data->current_iterate().delta_upper() * data->current_iterate().dual_bounds_u() -mu,
                    VecRealScalar(m, 0.)) +
            data->current_iterate().delta_upper() * data->current_iterate().delta_dual_bounds_u() +
            -1. * data->current_iterate().dual_bounds_u() *
                data->current_iterate().delta_primal_s()),
        1e-6);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
