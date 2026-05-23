//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#include "fatrop/graph/hessian.hpp"
#include "fatrop/graph/jacobian.hpp"
#include "fatrop/graph/nlp_graph.hpp"
#include "fatrop/graph/problem_info.hpp"
#include "fatrop/graph/problem_type.hpp"
#include "fatrop/ip_algorithm/ip_iterate.hxx"

using namespace fatrop;

// For graph problems offset_g_eq_slack == 0 and number_of_g_eq_slack ==
// number_of_eq_constraints, so the inequality-violation slice is the entire
// constraint-violation vector.
template <> const VecRealView IpIterate<GraphProblem>::constr_viol_ineq()
{
    return constr_viol().block(info_->number_of_g_eq_slack, info_->offset_g_eq_slack);
}

// explicit template instantiation
template class fatrop::IpIterate<GraphProblem>;
