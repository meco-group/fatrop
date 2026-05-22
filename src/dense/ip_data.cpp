#include "fatrop/dense/hessian.hpp"
#include "fatrop/dense/jacobian.hpp"
#include "fatrop/dense/nlp_dense.hpp"
#include "fatrop/dense/problem_info.hpp"
#include "fatrop/dense/type.hpp"
#include "fatrop/ip_algorithm/ip_data.hxx"
#include "fatrop/ip_algorithm/ip_iterate.hpp"

using namespace fatrop;
// explicit template instantiation
template class fatrop::IpData<DenseType>;
