//
// Copyright (c) Lander Vanroye, KU Leuven
//
#include "fatrop/ocp/hessian.hpp"
#include "fatrop/common/exception.hpp"
#include "fatrop/linear_algebra/linear_algebra.hpp"
#include "fatrop/ocp/dims.hpp"
#include "fatrop/ocp/problem_info.hpp"
using namespace fatrop;

Hessian<OcpType>::Hessian(const ProblemDims<OcpType> &dims)
{
    // reserve memory for the Jacobian matrices
    RSQrqt.reserve(dims.K);
    // The Lagrangian Hessian lives in tangent space (it multiplies delta_u, delta_x).
    for (Index k = 0; k < dims.K; ++k)
        RSQrqt.emplace_back(
            dims.number_of_tangent_states[k] + dims.number_of_tangent_controls[k] + 1,
            dims.number_of_tangent_states[k] + dims.number_of_tangent_controls[k]);
};
void Hessian<OcpType>::apply_on_right(const OcpInfo &info, const VecRealView &x, Scalar alpha,
                                      const VecRealView &y, VecRealView &out) const
{
    for (Index k = 0; k < info.dims.K; ++k)
    {
        // Tangent-space dimensions: RSQ rows/cols match the search-direction layout.
        Index nu = info.dims.number_of_tangent_controls[k];
        Index nx = info.dims.number_of_tangent_states[k];
        // x, y and out are tangent-sized vectors here.
        Index offset_ux = info.offsets_tangent_u[k];
        // apply out[offs:offs+nu+nx] =  RSQ @ x[offs:offs+nu+nx]
        gemv_t(nu + nx, nu + nx, 1.0, RSQrqt[k], 0, 0, x, offset_ux, alpha, y, offset_ux, out,
               offset_ux);
    }
};
void Hessian<OcpType>::get_rhs(const OcpInfo &info, VecRealView &out) const
{
    for (Index k = 0; k < info.dims.K; ++k)
    {
        Index nu = info.dims.number_of_tangent_controls[k];
        Index nx = info.dims.number_of_tangent_states[k];
        Index offset_ux = info.offsets_tangent_u[k];
        // the rhs is the last row of the RSQrqt[k] matrix
        rowex(nu + nx, 1.0, RSQrqt[k], nu + nx, 0, out, offset_ux);
    }
};
void Hessian<OcpType>::set_rhs(const OcpInfo &info, const VecRealView &in)
{
    for (Index k = 0; k < info.dims.K; ++k)
    {
        Index nu = info.dims.number_of_tangent_controls[k];
        Index nx = info.dims.number_of_tangent_states[k];
        Index offset_ux = info.offsets_tangent_u[k];
        // the rhs is the last row of the RSQrqt[k] matrix
        rowin(nu + nx, 1.0, in, offset_ux, RSQrqt[k], nu + nx, 0);
    }
};
void Hessian<OcpType>::set_zero()
{
    for (auto &RSQ : RSQrqt)
        gese(RSQ.m(), RSQ.n(), 0.0, RSQ, 0, 0);
}
// make printable
namespace fatrop
{

    std::ostream &operator<<(std::ostream &os, const Hessian<OcpType> &hess)
    {
        os << "Hessian<OcpType> object with horizon length " << hess.RSQrqt.size();
        for (int k = 0; k < hess.RSQrqt.size(); ++k)
        {
            os << "\n ----- Stage " << k << ": -----\n";
            os << "RSQrq:\n" << transpose(hess.RSQrqt[k]) << "\n";
        }
        return os;
    }
}