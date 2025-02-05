#include "fatrop/ip_algorithm/ip_linesearch.hxx"
#include "fatrop/ocp/type.hpp"
#include "fatrop/ocp/problem_info.hpp"
#include "fatrop/ocp/jacobian.hpp"
#include "fatrop/ocp/hessian.hpp"
#include "fatrop/ocp/pd_system_orig.hpp"
#include "fatrop/ocp/pd_solver_orig.hpp"
#include "fatrop/ocp/nlp_ocp.hpp"
using namespace fatrop;
template class fatrop::IpLinesearch<OcpType>;