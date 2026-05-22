//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//
#include "fatrop/dense/jacobian.hpp"
#include "fatrop/dense/dims.hpp"
#include "fatrop/dense/problem_info.hpp"
#include "fatrop/linear_algebra/linear_algebra.hpp"

#include <ostream>

using namespace fatrop;

Jacobian<DenseType>::Jacobian(const ProblemDims<DenseType> &dims)
    : Gg_eqt(dims.nx_tangent + 1, dims.ng), Gg_ineqt(dims.nx_tangent + 1, dims.ng_ineq)
{
}

void Jacobian<DenseType>::apply_on_right(const DenseInfo &info, const VecRealView &x, Scalar alpha,
                                         const VecRealView &y, VecRealView &out) const
{
    const Index nxt = info.dims.nx_tangent;
    const Index ng = info.dims.ng;
    const Index ng_ineq = info.dims.ng_ineq;
    out = alpha * y;
    // The path-eq block sits at the front of the equality vector, the slack-eq block
    // follows it at offset info.offset_g_eq_slack.
    if (ng > 0)
        gemv_t(nxt, ng, 1.0, Gg_eqt, 0, 0, x, 0, 1.0, out, 0, out, 0);
    if (ng_ineq > 0)
        gemv_t(nxt, ng_ineq, 1.0, Gg_ineqt, 0, 0, x, 0, 1.0, out, info.offset_g_eq_slack, out,
               info.offset_g_eq_slack);
}

void Jacobian<DenseType>::transpose_apply_on_right(const DenseInfo &info,
                                                   const VecRealView &mult_eq, Scalar alpha,
                                                   const VecRealView &y, VecRealView &out) const
{
    const Index nxt = info.dims.nx_tangent;
    const Index ng = info.dims.ng;
    const Index ng_ineq = info.dims.ng_ineq;
    out = alpha * y;
    if (ng > 0)
        gemv_n(nxt, ng, 1.0, Gg_eqt, 0, 0, mult_eq, 0, 1.0, out, 0, out, 0);
    if (ng_ineq > 0)
        gemv_n(nxt, ng_ineq, 1.0, Gg_ineqt, 0, 0, mult_eq, info.offset_g_eq_slack, 1.0, out, 0,
               out, 0);
}

void Jacobian<DenseType>::get_rhs(const DenseInfo &info, VecRealView &rhs) const
{
    const Index nxt = info.dims.nx_tangent;
    if (info.dims.ng > 0)
        rowex(info.dims.ng, 1.0, Gg_eqt, nxt, 0, rhs, 0);
    if (info.dims.ng_ineq > 0)
        rowex(info.dims.ng_ineq, 1.0, Gg_ineqt, nxt, 0, rhs, info.offset_g_eq_slack);
}

void Jacobian<DenseType>::set_rhs(const DenseInfo &info, const VecRealView &rhs)
{
    const Index nxt = info.dims.nx_tangent;
    if (info.dims.ng > 0)
        rowin(info.dims.ng, 1.0, rhs, 0, Gg_eqt, nxt, 0);
    if (info.dims.ng_ineq > 0)
        rowin(info.dims.ng_ineq, 1.0, rhs, info.offset_g_eq_slack, Gg_ineqt, nxt, 0);
}

namespace fatrop
{
    std::ostream &operator<<(std::ostream &os, const Jacobian<DenseType> &jac)
    {
        os << "Jacobian<DenseType>\n";
        os << "Gg_eq:\n" << transpose(jac.Gg_eqt) << "\n";
        os << "Gg_ineq:\n" << transpose(jac.Gg_ineqt) << "\n";
        return os;
    }
} // namespace fatrop
