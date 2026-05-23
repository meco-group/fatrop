//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#include "fatrop/graph/hessian.hpp"
#include "fatrop/graph/jacobian.hpp"
#include "fatrop/graph/problem_info.hpp"
#include "fatrop/graph/problem_type.hpp"
#include "fatrop/ip_algorithm/ip_resto_phase_min_cl1.hxx"

template class fatrop::IpRestoPhaseMinCl1<fatrop::GraphProblem>;
