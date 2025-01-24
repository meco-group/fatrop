#include "fatrop/ip_algorithm/ip_eq_mult_initializer.hxx"
#include "fatrop/ocp/hessian.hpp"
#include "fatrop/ocp/jacobian.hpp"
#include "fatrop/ocp/nlp_ocp.hpp"
#include "fatrop/ocp/pd_solver_orig.hpp"
#include "fatrop/ocp/pd_system_orig.hpp"
#include "fatrop/ocp/type.hpp"

using namespace fatrop;

template class IpEqMultInitializer<OcpType>; 