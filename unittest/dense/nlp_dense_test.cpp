//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//
#include "fatrop/dense/dense_abstract.hpp"
#include "fatrop/dense/hessian.hpp"
#include "fatrop/dense/jacobian.hpp"
#include "fatrop/dense/nlp_dense.hpp"
#include "fatrop/dense/problem_info.hpp"
#include "fatrop/dense/type.hpp"
#include "fatrop/linear_algebra/linear_algebra.hpp"
#include "fatrop/nlp/dims.hpp"
#include "dense_test_problem.hpp"
#include <gtest/gtest.h>
#include <memory>

using namespace fatrop;
using namespace fatrop::test;

class NlpDenseTest : public ::testing::Test
{
protected:
    std::shared_ptr<DenseTestProblem> problem = std::make_shared<DenseTestProblem>();
    NlpDense nlp{problem};
    ProblemInfo<DenseType> info{nlp.problem_dims()};
    VecRealAllocated primal_x{nlp.nlp_dims().number_of_variables};
    VecRealAllocated primal_s{nlp.nlp_dims().number_of_ineq_constraints};
    VecRealAllocated mult{nlp.nlp_dims().number_of_eq_constraints};

    void SetUp() override
    {
        nlp.get_initial_primal(info, primal_x);
        primal_s = 0.1;
        mult = 0.;
    }
};

TEST_F(NlpDenseTest, CheckProblemDimensions)
{
    EXPECT_EQ(nlp.problem_dims().nx, DenseTestProblemMath::nx);
    EXPECT_EQ(nlp.problem_dims().ng, DenseTestProblemMath::ng);
    EXPECT_EQ(nlp.problem_dims().ng_ineq, DenseTestProblemMath::ng_ineq);
    EXPECT_EQ(nlp.nlp_dims().number_of_variables, DenseTestProblemMath::nx);
    EXPECT_EQ(nlp.nlp_dims().number_of_tangent_variables, DenseTestProblemMath::nx);
    EXPECT_EQ(nlp.nlp_dims().number_of_eq_constraints,
              DenseTestProblemMath::ng + DenseTestProblemMath::ng_ineq);
    EXPECT_EQ(nlp.nlp_dims().number_of_ineq_constraints, DenseTestProblemMath::ng_ineq);
}

TEST_F(NlpDenseTest, TestJacobian)
{
    Jacobian<DenseType> jac(nlp.problem_dims());
    EXPECT_NO_THROW(nlp.eval_constr_jac(info, primal_x, primal_s, jac));
}

TEST_F(NlpDenseTest, TestLagrangianHessian)
{
    Hessian<DenseType> hess(nlp.problem_dims());
    Scalar objective_scale = 1.0;
    EXPECT_NO_THROW(nlp.eval_lag_hess(info, objective_scale, primal_x, primal_s, mult, hess));
}

TEST_F(NlpDenseTest, TestConstraintViolation)
{
    VecRealAllocated res(nlp.nlp_dims().number_of_eq_constraints);
    EXPECT_NO_THROW(nlp.eval_constraint_violation(info, primal_x, primal_s, res));
}

TEST_F(NlpDenseTest, TestObjectiveGradient)
{
    VecRealAllocated grad_x(nlp.nlp_dims().number_of_tangent_variables);
    VecRealAllocated grad_s(nlp.nlp_dims().number_of_ineq_constraints);
    Scalar objective_scale = 1.0;
    EXPECT_NO_THROW(
        nlp.eval_objective_gradient(info, objective_scale, primal_x, primal_s, grad_x, grad_s));
}

TEST_F(NlpDenseTest, TestObjective)
{
    Scalar val = 0.0;
    Scalar objective_scale = 1.0;
    EXPECT_NO_THROW(nlp.eval_objective(info, objective_scale, primal_x, primal_s, val));
    // Closed form: f(x_init = (0.1, 0.2, 0.3, 0.4)).
    Scalar expected = 0.5 * (1. * 0.01 + 2. * 0.04 + 3. * 0.09 + 4. * 0.16) +
                      (-1.) * 0.1 + (-2.) * 0.2 + (-3.) * 0.3 + (-4.) * 0.4;
    EXPECT_NEAR(val, expected, 1e-12);
}
