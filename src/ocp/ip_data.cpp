#include "fatrop/ip_algorithm/ip_iterate.hpp"
#include "fatrop/ip_algorithm/ip_data.hxx"
#include "fatrop/ocp/nlp_ocp.hpp"
#include "fatrop/ocp/type.hpp"
#include "fatrop/ocp/problem_info.hpp"
#include "fatrop/ocp/jacobian.hpp"
#include "fatrop/ocp/hessian.hpp"

using namespace fatrop;
// explicit template instantiation
template class fatrop::IpData<OcpType>; 