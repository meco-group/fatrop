#include "fatrop/dense/hessian.hpp"
#include "fatrop/dense/jacobian.hpp"
#include "fatrop/dense/problem_info.hpp"
#include "fatrop/dense/type.hpp"
#include "fatrop/ip_algorithm/ip_convergence_check.hxx"

using namespace fatrop;

template class fatrop::IpConvergenceCheck<DenseType>;
