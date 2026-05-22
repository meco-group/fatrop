//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//
#include "fatrop/dense/hessian.hpp"
#include "fatrop/dense/dims.hpp"
#include "fatrop/dense/problem_info.hpp"
#include "fatrop/linear_algebra/linear_algebra.hpp"

#include <ostream>

using namespace fatrop;

Hessian<DenseType>::Hessian(const ProblemDims<DenseType> &dims)
    : Hht(dims.nx_tangent + 1, dims.nx_tangent)
{
}

void Hessian<DenseType>::apply_on_right(const DenseInfo &info, const VecRealView &x, Scalar alpha,
                                        const VecRealView &y, VecRealView &out) const
{
    const Index nxt = info.dims.nx_tangent;
    gemv_t(nxt, nxt, 1.0, Hht, 0, 0, x, 0, alpha, y, 0, out, 0);
}

void Hessian<DenseType>::get_rhs(const DenseInfo &info, VecRealView &out) const
{
    const Index nxt = info.dims.nx_tangent;
    rowex(nxt, 1.0, Hht, nxt, 0, out, 0);
}

void Hessian<DenseType>::set_rhs(const DenseInfo &info, const VecRealView &in)
{
    const Index nxt = info.dims.nx_tangent;
    rowin(nxt, 1.0, in, 0, Hht, nxt, 0);
}

void Hessian<DenseType>::set_zero() { gese(Hht.m(), Hht.n(), 0.0, Hht, 0, 0); }

namespace fatrop
{
    std::ostream &operator<<(std::ostream &os, const Hessian<DenseType> &hess)
    {
        os << "Hessian<DenseType>\n";
        os << "Hh:\n" << transpose(hess.Hht) << "\n";
        return os;
    }
} // namespace fatrop
