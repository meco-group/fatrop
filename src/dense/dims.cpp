//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//
#include "fatrop/dense/dims.hpp"
#include "fatrop/common/exception.hpp"

using namespace fatrop;

ProblemDims<DenseType>::ProblemDims(Index nx, Index ng, Index ng_ineq)
    : nx(nx), nx_tangent(nx), ng(ng), ng_ineq(ng_ineq)
{
    check_problem_dimensions();
}

ProblemDims<DenseType>::ProblemDims(Index nx, Index nx_tangent, Index ng, Index ng_ineq)
    : nx(nx), nx_tangent(nx_tangent), ng(ng), ng_ineq(ng_ineq)
{
    check_problem_dimensions();
}

void ProblemDims<DenseType>::check_problem_dimensions() const
{
    fatrop_assert_msg(nx >= 0, "The number of primal variables must be non-negative.");
    fatrop_assert_msg(nx_tangent >= 0, "The tangent dimension must be non-negative.");
    fatrop_assert_msg(ng >= 0, "The number of equality constraints must be non-negative.");
    fatrop_assert_msg(ng_ineq >= 0, "The number of inequality constraints must be non-negative.");
    fatrop_assert_msg(ng <= nx_tangent,
                      "The number of equality constraints exceeds the tangent dimension.");
}
