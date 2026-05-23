//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#include "fatrop/graph/problem_info.hpp"

using namespace fatrop;

ProblemInfo<GraphProblem>::ProblemInfo(const ProblemDims<GraphProblem> &dims_in)
    : dims(dims_in), sparsity(dims_in.sparsity_ptr)
{
    number_of_primal_variables = dims.nx;
    number_of_tangent_variables = dims.nx_tan;

    number_of_slack_variables = dims.ng_ineq_total;
    // D_x layout: [primal-tangent damping | slack damping].
    offset_slack = number_of_tangent_variables;

    // No equality (path) constraints — only the slack-equality block is
    // present.
    number_of_g_eq_slack = dims.ng_ineq_total;
    number_of_eq_constraints = number_of_g_eq_slack;
    offset_g_eq_slack = 0;

    // PD system (original) layout: [primal-tangent | slack | mult | zl | zu].
    pd_orig_offset_primal = 0;
    pd_orig_offset_slack = pd_orig_offset_primal + number_of_tangent_variables;
    pd_orig_offset_mult = pd_orig_offset_slack + number_of_slack_variables;
    pd_orig_offset_zl = pd_orig_offset_mult + number_of_eq_constraints;
    pd_orig_offset_zu = pd_orig_offset_zl + number_of_slack_variables;

    // Graph problems have no restoration phase. The following fields are set
    // to consistent values mirroring the dense layout so the (untaken)
    // shared-code branches that touch them do not produce garbage.
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

    // Graph problems only have inequality-derived equality constraints (the
    // slack-eq rows); dual damping applies uniformly.
    constraint_allows_dual_damping = std::vector<bool>(number_of_eq_constraints, true);
}
