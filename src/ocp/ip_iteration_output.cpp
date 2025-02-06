#include "fatrop/ip_algorithm/ip_iteration_output.hxx"
#include "fatrop/ocp/type.hpp"
// todo avoid that i have to include this here
#include "fatrop/ocp/jacobian.hpp"
#include "fatrop/ocp/hessian.hpp"
#include "fatrop/ocp/problem_info.hpp"

template class fatrop::IpIterationOutput<fatrop::OcpType>;
