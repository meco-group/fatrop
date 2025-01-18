#include "fatrop/ip_algorithm/ip_iterate.hxx"
#include "fatrop/ocp/nlp_ocp.hpp"
#include "fatrop/ocp/type.hpp"
#include "fatrop/ocp/jacobian.hpp"
#include "fatrop/ocp/hessian.hpp"

using namespace fatrop;
// explicit template instantiation
template class IpIterate<OcpType>; 