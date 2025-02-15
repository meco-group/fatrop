#include "fatrop/ip_algorithm/ip_initializer_resto.hxx"
#include "fatrop/ocp/type.hpp"
#include "fatrop/ocp/problem_info.hpp"
#include "fatrop/ocp/hessian.hpp"
#include "fatrop/ocp/jacobian.hpp"

using namespace fatrop;
// explicit template instantiation
template class fatrop::IpInitializerResto<OcpType>;