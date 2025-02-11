#include "fatrop/ip_algorithm/ip_alg_builder.hxx"
#include "fatrop/ocp/aug_system_solver.hpp"
#include "fatrop/ocp/hessian.hpp"
#include "fatrop/ocp/jacobian.hpp"
#include "fatrop/ocp/nlp_ocp.hpp"
#include "fatrop/ocp/pd_solver_orig.hpp"
#include "fatrop/ocp/pd_system_orig.hpp"
#include "fatrop/ocp/pd_solver_resto.hpp"
#include "fatrop/ocp/pd_system_resto.hpp"
#include "fatrop/ocp/problem_info.hpp"
#include "fatrop/ocp/type.hpp"
#include "fatrop/ocp/aug_system_solver.hpp"
// #include "fatrop/ocp/problem_info.hpp"

namespace fatrop
{
    template class IpAlgBuilder<OcpType>;
} // namespace fatrop
