//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//
#include "fatrop/dense/problem_info.hpp"

using namespace fatrop;

ProblemInfo<DenseType>::ProblemInfo(const ProblemDims<DenseType> &dims) : dims(dims)
{
    number_of_primal_variables = dims.nx;
    number_of_tangent_variables = dims.nx_tangent;

    number_of_slack_variables = dims.ng_ineq;
    // D_x is laid out as [primal-tangent damping | slack damping] inside the
    // IP iterate, so the slack block starts after the tangent-space block.
    offset_slack = number_of_tangent_variables;

    number_of_g_eq_slack = dims.ng_ineq;
    number_of_eq_constraints = dims.ng + number_of_g_eq_slack;
    // Equality vector is laid out as [g_eq_path, g_eq_slack]; the slack block follows
    // the (ng-long) path block.
    offset_g_eq_slack = dims.ng;

    // PD system (original) layout: [primal-tangent | slack | mult | zl | zu].
    // The primal block stores the search direction (tangent-space), so its size
    // is number_of_tangent_variables — identical to number_of_primal_variables
    // for Euclidean problems.
    pd_orig_offset_primal = 0;
    pd_orig_offset_slack = pd_orig_offset_primal + number_of_tangent_variables;
    pd_orig_offset_mult = pd_orig_offset_slack + number_of_slack_variables;
    pd_orig_offset_zl = pd_orig_offset_mult + number_of_eq_constraints;
    pd_orig_offset_zu = pd_orig_offset_zl + number_of_slack_variables;

    // Restoration phase adds n and p slacks per equality constraint.
    number_of_slack_variables_resto = number_of_slack_variables + 2 * number_of_eq_constraints;
    pd_resto_offset_primal = 0;
    pd_resto_offset_slack = pd_resto_offset_primal + number_of_tangent_variables;
    pd_resto_offset_mult = pd_resto_offset_slack + number_of_slack_variables_resto;
    pd_resto_offset_zl = pd_resto_offset_mult + number_of_eq_constraints;
    pd_resto_offset_zu = pd_resto_offset_zl + number_of_slack_variables_resto;
    pd_resto_offset_zp = pd_resto_offset_zl + number_of_slack_variables;
    pd_resto_offset_zn = pd_resto_offset_zp + number_of_eq_constraints;

    offset_slack_p = offset_slack + number_of_slack_variables;
    offset_slack_n = offset_slack_p + number_of_eq_constraints;
    offset_s = 0;
    offset_p = number_of_slack_variables;
    offset_n = offset_p + number_of_eq_constraints;

    // No dynamics rows to exclude, so all equality constraints allow dual damping.
    constraint_allows_dual_damping = std::vector<bool>(number_of_eq_constraints, true);
}
