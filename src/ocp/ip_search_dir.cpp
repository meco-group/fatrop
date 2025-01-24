#include "fatrop/ip_algorithm/ip_search_dir.hxx"
#include "fatrop/ocp/nlp_ocp.hpp"
#include "fatrop/ocp/type.hpp"
#include "fatrop/ocp/jacobian.hpp"
#include "fatrop/ocp/hessian.hpp"
#include "fatrop/ocp/pd_solver_orig.hpp"
#include "fatrop/ocp/pd_system_orig.hpp"

using namespace fatrop;
// explicit template instantiation
template class IpSearchDirImpl<OcpType>; 