#include "fatrop/ip_algorithm/ip_convergence_check_resto.hxx"
#include "fatrop/ip_algorithm/ip_convergence_check.hxx"
#include "fatrop/ocp/type.hpp"
#include "fatrop/ocp/problem_info.hpp"
#include "fatrop/ocp/jacobian.hpp"
#include "fatrop/ocp/hessian.hpp"
using namespace fatrop;

// explicit template instantiation
template class fatrop::IpConvergenceCheckResto<OcpType>;