//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#include "fatrop/graph/jacobian.hpp"

#include "fatrop/graph/problem_info.hpp"
#include "fatrop/linear_algebra/blasfeo_operations.hpp"

#include <ostream>

using namespace fatrop;

namespace
{
// Build per-block transposed-Jacobian matrices with one trailing row reserved
// for the right-hand side (block-local g_ineq_k value). For blocks with no
// inequality constraints we still allocate a 1-column matrix to keep
// per-block bookkeeping uniform — the column is never accessed.
std::vector<MatRealAllocated> make_block_ineq_jac(const ProblemDims<GraphProblem> &dims)
{
    std::vector<MatRealAllocated> out;
    out.reserve(dims.num_blocks);
    for (Index k = 0; k < dims.num_blocks; ++k)
    {
        const Index ngk = std::max<Index>(dims.ng_ineq[k], 1);
        // Tangent-space Jacobian columns + RHS row.
        out.emplace_back(dims.block_sizes_tan[k] + 1, ngk);
    }
    return out;
}
} // namespace

Jacobian<GraphProblem>::Jacobian(const ProblemDims<GraphProblem> &dims)
    : Gg_ineqt(make_block_ineq_jac(dims))
{
    // Zero the trailing RHS row of every block; the column data is filled in
    // by NlpGraph::eval_constr_jac.
    for (auto &M : Gg_ineqt)
        gese(M.m(), M.n(), 0.0, M, 0, 0);
}

void Jacobian<GraphProblem>::apply_on_right(const GraphInfo &info, const VecRealView &x,
                                            Scalar alpha, const VecRealView &y,
                                            VecRealView &out) const
{
    const Index nx = info.dims.nx;
    const Index ng_total = info.dims.ng_ineq_total;
    (void)nx;
    // out = alpha * y, then add block-local A_k * x_k for each block.
    if (&out != &y || alpha != 1.0)
    {
        out = alpha * y;
    }
    for (Index k = 0; k < info.dims.num_blocks; ++k)
    {
        const Index nk = info.dims.block_sizes_tan[k];
        const Index ngk = info.dims.ng_ineq[k];
        if (ngk == 0)
            continue;
        const Index x_off = info.dims.block_offsets_tan[k];
        const Index ineq_off = info.offset_g_eq_slack + info.dims.ng_ineq_offsets[k];
        gemv_t(nk, ngk, 1.0, Gg_ineqt[k], 0, 0, x, x_off, 1.0, out, ineq_off, out, ineq_off);
    }
    (void)ng_total;
}

void Jacobian<GraphProblem>::transpose_apply_on_right(const GraphInfo &info,
                                                      const VecRealView &mult_eq, Scalar alpha,
                                                      const VecRealView &y, VecRealView &out) const
{
    if (&out != &y || alpha != 1.0)
    {
        out = alpha * y;
    }
    for (Index k = 0; k < info.dims.num_blocks; ++k)
    {
        const Index nk = info.dims.block_sizes_tan[k];
        const Index ngk = info.dims.ng_ineq[k];
        if (ngk == 0)
            continue;
        const Index x_off = info.dims.block_offsets_tan[k];
        const Index ineq_off = info.offset_g_eq_slack + info.dims.ng_ineq_offsets[k];
        gemv_n(nk, ngk, 1.0, Gg_ineqt[k], 0, 0, mult_eq, ineq_off, 1.0, out, x_off, out, x_off);
    }
}

void Jacobian<GraphProblem>::get_rhs(const GraphInfo &info, VecRealView &rhs) const
{
    for (Index k = 0; k < info.dims.num_blocks; ++k)
    {
        const Index ngk = info.dims.ng_ineq[k];
        if (ngk == 0)
            continue;
        const Index nk = info.dims.block_sizes_tan[k];
        const Index ineq_off = info.offset_g_eq_slack + info.dims.ng_ineq_offsets[k];
        rowex(ngk, 1.0, Gg_ineqt[k], nk, 0, rhs, ineq_off);
    }
}

void Jacobian<GraphProblem>::set_rhs(const GraphInfo &info, const VecRealView &rhs)
{
    for (Index k = 0; k < info.dims.num_blocks; ++k)
    {
        const Index ngk = info.dims.ng_ineq[k];
        if (ngk == 0)
            continue;
        const Index nk = info.dims.block_sizes_tan[k];
        const Index ineq_off = info.offset_g_eq_slack + info.dims.ng_ineq_offsets[k];
        rowin(ngk, 1.0, rhs, ineq_off, Gg_ineqt[k], nk, 0);
    }
}

namespace fatrop
{
    std::ostream &operator<<(std::ostream &os, const Jacobian<GraphProblem> &jac)
    {
        os << "Jacobian<GraphProblem> (" << jac.Gg_ineqt.size() << " block-local rows)\n";
        return os;
    }
} // namespace fatrop
