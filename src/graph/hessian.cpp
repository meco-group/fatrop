//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#include "fatrop/graph/hessian.hpp"

#include "fatrop/graph/problem_info.hpp"
#include "fatrop/linear_algebra/blasfeo_operations.hpp"

#include <ostream>

using namespace fatrop;

Hessian<GraphProblem>::Hessian(const ProblemDims<GraphProblem> &dims)
    : H(dims.sparsity()), grad(dims.nx_tan)
{
    H.set_zero();
    grad = 0.0;
}

void Hessian<GraphProblem>::apply_on_right(const GraphInfo & /*info*/, const VecRealView &x,
                                            Scalar alpha, const VecRealView &y,
                                            VecRealView &out) const
{
    H.apply_on_right(x, alpha, y, out);
}

void Hessian<GraphProblem>::get_rhs(const GraphInfo &info, VecRealView &out) const
{
    veccp(info.dims.nx_tan, grad, 0, out, 0);
}

void Hessian<GraphProblem>::set_rhs(const GraphInfo &info, const VecRealView &in)
{
    veccp(info.dims.nx_tan, in, 0, grad, 0);
}

void Hessian<GraphProblem>::set_zero()
{
    H.set_zero();
    grad = 0.0;
}

namespace fatrop
{
    std::ostream &operator<<(std::ostream &os, const Hessian<GraphProblem> &hess)
    {
        os << "Hessian<GraphProblem> (block-sparse, " << hess.H.sparsity().num_blocks()
           << " blocks)\n";
        os << "grad: ";
        for (Index i = 0; i < hess.grad.m(); ++i)
            os << hess.grad(i) << ' ';
        os << "\n";
        return os;
    }
} // namespace fatrop
