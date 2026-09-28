// Copyright (c) 2026 Lander Vanroye, KU Leuven
#ifndef FATROP_GRAPH_BFGS_MATRIX_OPS_HPP
#define FATROP_GRAPH_BFGS_MATRIX_OPS_HPP

#include "fatrop/graph/hessian.hpp"
#include "fatrop/nlp/bfgs_matrix_ops.hpp"
#include <algorithm>
#include <cmath>
#include <utility>

namespace fatrop
{
    template <> struct BfgsMatrixOps<GraphProblem>
    {
        template <typename F> static void for_each_block(const Hessian<GraphProblem> &h, F &&f)
        { f(0, 0, h.H.sparsity().total_size()); }

        template <typename F> static void for_each_lower(
            const Hessian<GraphProblem> &h, Index, Index, F &&f)
        {
            const auto &sp = h.H.sparsity();
            for (Index j = 0; j < sp.num_blocks(); ++j)
                for (Index i : sp.column_pattern(j))
                {
                    auto entry = h.H.block(i, j);
                    for (Index c = 0; c < entry.n(); ++c)
                        for (Index r = i == j ? c : 0; r < entry.m(); ++r)
                            f(sp.block_offset(i) + r, sp.block_offset(j) + c, entry(r, c));
                }
        }

        static void set(Hessian<GraphProblem> &h, Index, Index,
                        Index i, Index j, Scalar value)
        {
            const auto &sp = h.H.sparsity();
            Index bi = block_of(sp, i), bj = block_of(sp, j);
            if (bi < bj) { std::swap(bi, bj); std::swap(i, j); }
            auto entry = h.H.block(bi, bj);
            entry(i - sp.block_offset(bi), j - sp.block_offset(bj)) = value;
            if (bi == bj)
                entry(j - sp.block_offset(bj), i - sp.block_offset(bi)) = value;
        }

        static void stabilize(Hessian<GraphProblem> &h, VecRealView &row_sum)
        {
            // Sparse projection can lose positive definiteness. Strict diagonal
            // dominance gives the graph Cholesky solver an SPD approximation.
            row_sum = 0.;
            for_each_lower(h, 0, 0, [&](Index i, Index j, Scalar value) {
                if (i != j)
                {
                    const Scalar magnitude = std::abs(value);
                    row_sum(i) += magnitude;
                    row_sum(j) += magnitude;
                }
            });
            for (Index i = 0; i < row_sum.m(); ++i)
            {
                const Scalar minimum = row_sum(i) + 1e-8;
                if (diagonal(h, i) < minimum) set(h, 0, 0, i, i, minimum);
            }
        }

        static void copy(const Hessian<GraphProblem> &source, Hessian<GraphProblem> &target)
        {
            const auto &sp = source.H.sparsity();
            for (Index j = 0; j < sp.num_blocks(); ++j)
                for (Index i : sp.column_pattern(j))
                {
                    auto from = source.H.block(i, j);
                    auto to = target.H.block(i, j);
                    for (Index c = 0; c < from.n(); ++c)
                        for (Index r = 0; r < from.m(); ++r) to(r, c) = from(r, c);
                }
            target.grad = source.grad;
        }

    private:
        static Index block_of(const BlockSparsityPattern &sp, Index coordinate)
        {
            return Index(std::upper_bound(sp.block_offsets().begin(), sp.block_offsets().end(),
                                          coordinate) - sp.block_offsets().begin()) - 1;
        }
        static Scalar diagonal(const Hessian<GraphProblem> &h, Index coordinate)
        {
            const auto &sp = h.H.sparsity();
            const Index block = block_of(sp, coordinate);
            const Index local = coordinate - sp.block_offset(block);
            return h.H.block(block, block)(local, local);
        }
    };
}
#endif
