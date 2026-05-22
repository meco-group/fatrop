#include "fatrop/dense/hessian.hpp"
#include "fatrop/dense/jacobian.hpp"
#include "fatrop/dense/problem_info.hpp"
#include "fatrop/dense/type.hpp"
#include "fatrop/ip_algorithm/ip_convergence_check.hxx"
#include "fatrop/ip_algorithm/ip_convergence_check_resto.hxx"

using namespace fatrop;

template class fatrop::IpConvergenceCheckResto<DenseType>;
