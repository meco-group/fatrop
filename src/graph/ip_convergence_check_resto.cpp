//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#include "fatrop/graph/hessian.hpp"
#include "fatrop/graph/jacobian.hpp"
#include "fatrop/graph/problem_info.hpp"
#include "fatrop/graph/problem_type.hpp"
#include "fatrop/ip_algorithm/ip_convergence_check.hxx"
#include "fatrop/ip_algorithm/ip_convergence_check_resto.hxx"

using namespace fatrop;

template class fatrop::IpConvergenceCheckResto<GraphProblem>;
