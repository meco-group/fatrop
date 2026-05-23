//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#include "fatrop/graph/hessian.hpp"
#include "fatrop/graph/jacobian.hpp"
#include "fatrop/graph/nlp_graph.hpp"
#include "fatrop/graph/problem_info.hpp"
#include "fatrop/graph/problem_type.hpp"
#include "fatrop/ip_algorithm/ip_data.hxx"
#include "fatrop/ip_algorithm/ip_iterate.hpp"

using namespace fatrop;

// explicit template instantiation
template struct fatrop::IpData<GraphProblem>;
