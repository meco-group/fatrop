#include "fatrop/ocp/jacobian.hpp"
#include "fatrop/common/exception.hpp"
using namespace fatrop;

Jacobian<OcpType>::Jacobian(const OcpDims &dims)
{
    // reserve memory for the Jacobian matrices
    jac_dyn.reserve(dims.K-1);
    jac_eq.reserve(dims.K);
    jac_ineq.reserve(dims.K);
    // allocate memory for the Jacobian matrices
    for (int k = 0; k < dims.K-1; ++k)
        jac_dyn.emplace_back(dims.number_of_states[k] + dims.number_of_constrols[k] + 1, dims.number_of_states[k+1]);
    for (int k = 0; k < dims.K; ++k)
    {
        fatrop_assert(dims.number_of_eq_constraints[k] <= dims.number_of_constrols[k] + dims.number_of_states[k] && "Number of equality constraints cannot exceed number of variables at a stage.");
        jac_eq.emplace_back(dims.number_of_states[k] + dims.number_of_constrols[k] + 1, dims.number_of_eq_constraints[k]);
    }
    for (int k = 0; k < dims.K; ++k)
        jac_ineq.emplace_back(dims.number_of_states[k] + dims.number_of_constrols[k] + 1, dims.number_of_ineq_constraints[k]);
};