#include "fatrop/linear_algebra/linear_algebra.hpp"
#include "fatrop/nlp/dims.hpp"
#include "fatrop/ocp/hessian.hpp"
#include "fatrop/ocp/jacobian.hpp"
#include "fatrop/ocp/nlp_ocp.hpp"
#include "fatrop/ocp/ocp_abstract.hpp"
#include "ocp_test_probem.hpp"
#include <gtest/gtest.h>
#include <memory>

using namespace fatrop;
using namespace fatrop::test;

class NlpOcpTest : public ::testing::Test
{
    typedef std::shared_ptr<OcpTestProblem> OcpTestSp;

protected:
    OcpTestSp ocp = std::make_shared<OcpTestProblem>();
    NlpOcp nlp_ocp = NlpOcp(ocp);
    ProblemInfo<OcpType> ocp_info = ProblemInfo<OcpType>(nlp_ocp.problem_dims());
    VecRealAllocated primal_x = VecRealAllocated(nlp_ocp.nlp_dims().number_of_variables);
    VecRealAllocated primal_s = VecRealAllocated(nlp_ocp.nlp_dims().number_of_ineq_constraints);
    VecRealAllocated mult = VecRealAllocated(nlp_ocp.nlp_dims().number_of_eq_constraints);
    void SetUp() override {}
};

TEST_F(NlpOcpTest, CheckProblemDimensions)
{
    EXPECT_EQ(nlp_ocp.problem_dims().K, 100);
    EXPECT_EQ(nlp_ocp.problem_dims().number_of_states[0], 4);
    EXPECT_EQ(nlp_ocp.problem_dims().number_of_controls[0], 2);
    EXPECT_EQ(nlp_ocp.problem_dims().number_of_controls[99], 0);
    EXPECT_EQ(nlp_ocp.problem_dims().number_of_eq_constraints[0], 4);
    EXPECT_EQ(nlp_ocp.problem_dims().number_of_eq_constraints[99], 4);
    EXPECT_EQ(nlp_ocp.problem_dims().number_of_ineq_constraints[0], 2);
    EXPECT_EQ(nlp_ocp.problem_dims().number_of_ineq_constraints[99], 0);
}

TEST_F(NlpOcpTest, TestJacobian)
{
    Jacobian<OcpType> jac(nlp_ocp.problem_dims());
    EXPECT_NO_THROW(nlp_ocp.eval_constr_jac(ocp_info, primal_x, primal_s, jac));
}

TEST_F(NlpOcpTest, TestLagrangianHessian)
{
    Hessian<OcpType> hess(nlp_ocp.problem_dims());
    Scalar objective_scale = 1.0;
    EXPECT_NO_THROW(
        nlp_ocp.eval_lag_hess(ocp_info, objective_scale, primal_x, primal_s, mult, hess));
    // Add more specific checks here if needed
}

TEST_F(NlpOcpTest, TestConstraintViolation)
{
    VecRealAllocated res(nlp_ocp.nlp_dims().number_of_eq_constraints +
                         nlp_ocp.nlp_dims().number_of_ineq_constraints);
    EXPECT_NO_THROW(nlp_ocp.eval_constraint_violation(ocp_info, primal_x, primal_s, res));
    // Add more specific checks here if needed
}

TEST_F(NlpOcpTest, TestObjectiveGradient)
{
    VecRealAllocated grad_x(nlp_ocp.nlp_dims().number_of_variables);
    VecRealAllocated grad_s(nlp_ocp.nlp_dims().number_of_ineq_constraints);
    Scalar objective_scale = 1.0;
    EXPECT_NO_THROW(
        nlp_ocp.eval_objective_gradient(ocp_info, objective_scale, primal_x, grad_x, grad_s));
    // Add more specific checks here if needed
}

TEST_F(NlpOcpTest, TestObjective)
{
    Scalar objective_value = 0.0;
    Scalar objective_scale = 1.0;
    EXPECT_NO_THROW(
        nlp_ocp.eval_objective(ocp_info, objective_scale, primal_x, primal_s, objective_value));
    // Add more specific checks here if needed
}
