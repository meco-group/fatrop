#include "fatrop/ip_algorithm/ip_linesearch.hxx"
#include "fatrop/ocp/type.hpp"
#include "fatrop/ocp/problem_info.hpp"
#include "fatrop/ocp/jacobian.hpp"
#include "fatrop/ocp/hessian.hpp"
#include "fatrop/ocp/pd_system_orig.hpp"
#include "fatrop/ocp/pd_solver_orig.hpp"
#include "fatrop/ocp/pd_system_resto.hpp"
#include "fatrop/ocp/pd_solver_resto.hpp"
#include "fatrop/ocp/nlp_ocp.hpp"

using namespace fatrop;

// Template instantiation for the original problem type
template class fatrop::IpLinesearch<PdSolverOrig<OcpType>, OcpType>;

// Template instantiation for the restoration problem type
template class fatrop::IpLinesearch<PdSolverResto<OcpType>, OcpType>;
