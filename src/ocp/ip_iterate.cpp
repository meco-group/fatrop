#include "fatrop/ip_algorithm/ip_iterate.hxx"
#include "fatrop/ocp/hessian.hpp"
#include "fatrop/ocp/jacobian.hpp"
#include "fatrop/ocp/nlp_ocp.hpp"
#include "fatrop/ocp/type.hpp"

using namespace fatrop;
// specializations
template <> const VecRealView IpIterate<OcpType>::constr_viol_ineq() { return constr_viol().block(info_.number_of_g_eq_slack, info_.offset_g_eq_slack); }
// explicit template instantiation
template class fatrop::IpIterate<OcpType>;