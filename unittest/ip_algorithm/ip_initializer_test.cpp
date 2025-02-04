#include "../ocp/ocp_test_probem.hpp"
#include "fatrop/ip_algorithm/ip_data.hpp"
#include "fatrop/ip_algorithm/ip_eq_mult_initializer.hpp"
#include "fatrop/ip_algorithm/ip_initializer.hpp"
#include "fatrop/ocp/aug_system_solver.hpp"
#include "fatrop/ocp/hessian.hpp"
#include "fatrop/ocp/jacobian.hpp"
#include "fatrop/ocp/nlp_ocp.hpp"
#include "fatrop/ocp/pd_solver_orig.hpp"
#include "fatrop/ocp/type.hpp"
#include <gtest/gtest.h>

namespace fatrop
{
    namespace test
    {

        class IpInitializerTest : public ::testing::Test
        {
        protected:
            IpInitializerTest()
                : problem(std::make_shared<OcpTestProblem>()),
                  nlp(std::make_shared<NlpOcp>(problem)), info(nlp->problem_dims()),
                  ipdata(std::make_shared<IpData<OcpType>>(nlp)),
                  D_x(info.number_of_primal_variables + info.number_of_slack_variables),
                  D_eq(info.number_of_eq_constraints),
                  aug_solver(std::make_shared<AugSystemSolver<OcpType>>(info)),
                  linear_solver(std::make_shared<PdSolverOrig<OcpType>>(info, aug_solver)),
                  eq_mult_initializer(
                      std::make_shared<IpEqMultInitializer<OcpType>>(ipdata, linear_solver)),
                  initializer(std::make_shared<IpInitializer<OcpType>>(ipdata, eq_mult_initializer))
            {
                ipdata->current_iterate().set_dual_bounds_l(
                    VecRealScalar(ipdata->current_iterate().dual_bounds_l().m(), 1));
                ipdata->current_iterate().set_dual_bounds_u(
                    VecRealScalar(ipdata->current_iterate().dual_bounds_u().m(), 1));
            }

            std::shared_ptr<OcpTestProblem> problem;
            std::shared_ptr<NlpOcp> nlp;
            ProblemInfo<OcpType> info;
            std::shared_ptr<IpData<OcpType>> ipdata;
            VecRealAllocated D_x, D_eq;
            std::shared_ptr<AugSystemSolver<OcpType>> aug_solver;
            std::shared_ptr<PdSolverOrig<OcpType>> linear_solver;
            std::shared_ptr<IpEqMultInitializer<OcpType>> eq_mult_initializer;
            std::shared_ptr<IpInitializer<OcpType>> initializer;
        };

        TEST_F(IpInitializerTest, EqMultInitializerTest)
        {
            // Check that the equality multiplier initialization doesn't throw an exception
            EXPECT_NO_THROW({ initializer->initialize(); })
                << "Equality multiplier initialization should not throw an exception";
        }

        TEST_F(IpInitializerTest, BoundPushTestLow)
        {

            auto &current_iterate = ipdata->current_iterate();
            auto &primal_s = current_iterate.primal_s();
            auto &primal_x = current_iterate.primal_x();
            current_iterate.set_primal_x(VecRealScalar(primal_x.m(), -100.));
            Scalar norm_l1_cv = norm_l1(current_iterate.constr_viol_ineq());

            // Initialize the problem
            initializer->initialize();

            const auto &lower_bounds = current_iterate.lower_bounds();
            const auto &upper_bounds = current_iterate.upper_bounds();

            for (int i = 0; i < primal_s.m(); ++i)
            {
                if (current_iterate.lower_bounded()[i])
                {
                    EXPECT_GT(primal_s(i), lower_bounds(i))
                        << "Slack should be pushed away from lower bound at index " << i;
                }
                if (current_iterate.upper_bounded()[i])
                {
                    EXPECT_LT(primal_s(i), upper_bounds(i))
                        << "Slack should be pushed away from upper bound at index " << i;
                }
            }
            EXPECT_LT(norm_l1(current_iterate.constr_viol_ineq()), norm_l1_cv)
                << "Norm of inequality constraint violation should decrease";
        }

        TEST_F(IpInitializerTest, BoundPushTestHigh)
        {

            auto &current_iterate = ipdata->current_iterate();
            auto &primal_s = current_iterate.primal_s();
            auto &primal_x = current_iterate.primal_x();
            current_iterate.set_primal_x(VecRealScalar(primal_x.m(), 100.));
            Scalar norm_l1_cv = norm_l1(current_iterate.constr_viol_ineq());

            // Initialize the problem
            initializer->initialize();

            const auto &lower_bounds = current_iterate.lower_bounds();
            const auto &upper_bounds = current_iterate.upper_bounds();

            for (int i = 0; i < primal_s.m(); ++i)
            {
                if (current_iterate.lower_bounded()[i])
                {
                    EXPECT_GT(primal_s(i), lower_bounds(i))
                        << "Slack should be pushed away from lower bound at index " << i;
                }
                if (current_iterate.upper_bounded()[i])
                {
                    EXPECT_LT(primal_s(i), upper_bounds(i))
                        << "Slack should be pushed away from upper bound at index " << i;
                }
            }
            EXPECT_LT(norm_l1(current_iterate.constr_viol_ineq()), norm_l1_cv)
                << "Norm of inequality constraint violation should decrease";
        }

    } // namespace fatrop::test
} // namespace fatrop::test

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
