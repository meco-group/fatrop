//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#include "fatrop/graph/block_pd_matrix.hpp"

#include "fatrop/common/exception.hpp"
#include "fatrop/linear_algebra/blasfeo_operations.hpp"

using namespace fatrop;

BlockPdMatrix::BlockPdMatrix(const BlockSparsityPattern &sp)
    : sp_(sp)
{
    const Index N = sp.num_blocks();
    // Reserve exactly so push_back never reallocates.
    blocks_.reserve(sp.num_stored_blocks());
    for (Index j = 0; j < N; ++j)
        for (Index i : sp.column_pattern(j))
            blocks_.emplace_back(sp.block_size(i), sp.block_size(j));
}

bool BlockPdMatrix::has_block(Index i, Index j) const
{
    if (i < j)
        std::swap(i, j);
    return sp_.block_id(i, j) >= 0;
}

MatRealView BlockPdMatrix::block(Index i, Index j)
{
    const Index e = sp_.block_id(i, j);
    fatrop_dbg_assert(e >= 0 && "block not in sparsity pattern");
    return blocks_[e];
}

const MatRealView BlockPdMatrix::block(Index i, Index j) const
{
    const Index e = sp_.block_id(i, j);
    fatrop_dbg_assert(e >= 0 && "block not in sparsity pattern");
    return blocks_[e];
}

void BlockPdMatrix::set_zero()
{
    for (auto &b : blocks_)
        b = 0.0;
}

void BlockPdMatrix::apply_on_right(const VecRealView &x, Scalar alpha, const VecRealView &y,
                                   VecRealView &out) const
{
    // No heap allocation in this function.
    const Index N = sp_.num_blocks();
    const Index m = sp_.total_size();
    // out = alpha * y. Handle aliasing safely.
    if (&out == &y)
    {
        if (alpha != 1.0)
            vecsc(m, alpha, out, 0);
    }
    else
    {
        veccpsc(m, alpha, y, 0, out, 0);
    }

    // Diagonal block contribution. We deliberately read only the lower
    // triangle of every stored diagonal block — the upper triangle is treated
    // as symmetric and not touched, so the caller does not need to mirror.
    for (Index k = 0; k < N; ++k)
    {
        const Index nk = sp_.block_size(k);
        const Index off = sp_.block_offset(k);
        if (nk == 0)
            continue;
        const MatRealView Mkk = blocks_[sp_.column_offset(k)];
        for (Index ii = 0; ii < nk; ++ii)
        {
            out(off + ii) += Mkk(ii, ii) * x(off + ii);
            for (Index jj = 0; jj < ii; ++jj)
            {
                const Scalar mij = Mkk(ii, jj);
                out(off + ii) += mij * x(off + jj);
                out(off + jj) += mij * x(off + ii);
            }
        }
    }

    // Off-diagonal blocks (i > j) contribute symmetrically to rows i and j.
    for (Index j = 0; j < N; ++j)
    {
        const Index nj = sp_.block_size(j);
        const Index off_j = sp_.block_offset(j);
        const auto &col = sp_.column_pattern(j);
        // col[0] is the diagonal block, handled above.
        for (size_t a = 1; a < col.size(); ++a)
        {
            const Index i = col[a];
            const Index ni = sp_.block_size(i);
            const Index off_i = sp_.block_offset(i);
            const MatRealView Mij = blocks_[sp_.column_offset(j) + a];
            gemv_n(ni, nj, 1.0, Mij, 0, 0, x, off_j, 1.0, out, off_i, out, off_i);
            gemv_t(ni, nj, 1.0, Mij, 0, 0, x, off_i, 1.0, out, off_j, out, off_j);
        }
    }
}
