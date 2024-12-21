#include "fatrop/ocp/hessian.hpp"
#include "fatrop/common/exception.hpp"
using namespace fatrop;

Hessian<OcpType>::Hessian(const OcpDims &dims)
{
    // reserve memory for the Jacobian matrices
    hess.reserve(dims.K);
    // allocate memory for the Jacobian matrices
    for (Index k = 0; k < dims.K; ++k)
        hess.emplace_back(dims.number_of_states[k] + dims.number_of_constrols[k] + 1, dims.number_of_states[k] + dims.number_of_constrols[k]);
};