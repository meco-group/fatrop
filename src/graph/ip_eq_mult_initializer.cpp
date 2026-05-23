//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#include "fatrop/graph/hessian.hpp"
#include "fatrop/graph/jacobian.hpp"
#include "fatrop/graph/nlp_graph.hpp"
#include "fatrop/graph/pd_solver_orig.hpp"
#include "fatrop/graph/pd_system_orig.hpp"
#include "fatrop/graph/problem_type.hpp"
#include "fatrop/ip_algorithm/ip_eq_mult_initializer.hxx"

using namespace fatrop;

template class fatrop::IpEqMultInitializer<GraphProblem>;
