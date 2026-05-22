#include "fatrop/dense/hessian.hpp"
#include "fatrop/dense/jacobian.hpp"
#include "fatrop/dense/problem_info.hpp"
#include "fatrop/dense/type.hpp"
#include "fatrop/ip_algorithm/ip_mu_update.hxx"

template class fatrop::IpMonotoneMuUpdate<fatrop::DenseType>;
