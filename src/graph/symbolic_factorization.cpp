//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#include "fatrop/graph/symbolic_factorization.hpp"

#include "fatrop/graph/block_sparsity.hpp"
#include "fatrop/graph/elimination_order.hpp"

#include <algorithm>

using namespace fatrop;

BlockSymbolicFactorization::BlockSymbolicFactorization(const BlockSparsityPattern &sp,
                                                       const BlockEliminationOrder &order)
    : num_blocks_(sp.num_blocks()), parent_(sp.num_blocks(), -1),
      lower_pattern_(sp.num_blocks())
{
    const Index N = sp.num_blocks();
    const auto &inv = order.inv_perm();

    // Permuted adjacency in elimination-order coordinates.
    //   A_lower[j] = rows i > j with original (i, j) nonzero (column view).
    //   A_row[k]   = columns j < k with original (k, j) nonzero (row view).
    std::vector<std::vector<Index>> A_lower(N);
    std::vector<std::vector<Index>> A_row(N);
    for (Index v = 0; v < N; ++v)
    {
        const Index jp = inv[v];
        for (Index u : sp.neighbors(v))
        {
            const Index ip = inv[u];
            if (ip > jp)
            {
                A_lower[jp].push_back(ip);
                A_row[ip].push_back(jp);
            }
        }
    }

    // Build the elimination tree (Liu's algorithm — Davis "Direct Methods for
    // Sparse Linear Systems" §4.1). Process rows k = 0..N-1; for each
    // original nonzero (k, j) with j < k, walk j up the partial tree being
    // built (path-compressing the @c ancestor links along the way) until we
    // hit either a node whose ancestor link is unset or one that already
    // points to k.
    std::vector<Index> ancestor(N, -1);
    for (Index k = 0; k < N; ++k)
    {
        for (Index j : A_row[k])
        {
            Index r = j;
            while (ancestor[r] != -1 && ancestor[r] != k)
            {
                const Index t = ancestor[r];
                ancestor[r] = k;
                r = t;
            }
            if (ancestor[r] == -1)
            {
                ancestor[r] = k;
                parent_[r] = k;
            }
        }
    }

    // Build children list (children of node k are nodes c < k with parent c = k).
    std::vector<std::vector<Index>> children(N);
    for (Index c = 0; c < N; ++c)
    {
        if (parent_[c] != -1)
            children[parent_[c]].push_back(c);
    }

    // Compute column patterns:
    //   pattern(L, :, j) = original A_lower(:, j) ∪ ⋃_{c child of j} (pattern(L, :, c) ∩ rows > j)
    // Process columns in order (children are < parent).
    std::vector<Index> mark(N, -1);
    for (Index j = 0; j < N; ++j)
    {
        auto &col = lower_pattern_[j];
        mark[j] = j;
        auto add = [&](Index i)
        {
            if (mark[i] != j)
            {
                mark[i] = j;
                col.push_back(i);
            }
        };
        for (Index i : A_lower[j])
            add(i);
        // Inherited rows are > j except the parent slot j itself, which the
        // mark already excludes.
        for (Index c : children[j])
            for (Index i : lower_pattern_[c])
                add(i);
        std::sort(col.begin(), col.end());
    }

    // Fundamental-supernode partition. Walk columns left-to-right; column k
    // joins the previous column's supernode iff parent(k-1) = k AND the
    // patterns line up so that |lp(k)| = |lp(k-1)| - 1 (the parent condition
    // forces lp(k) ⊇ lp(k-1) \ {k}, so the count check seals equality).
    std::vector<Index> fundamental{0};
    for (Index k = 1; k < N; ++k)
    {
        const Index prev = k - 1;
        const bool can_merge =
            parent_[prev] == k &&
            static_cast<Index>(lower_pattern_[k].size()) + 1 ==
                static_cast<Index>(lower_pattern_[prev].size());
        if (!can_merge)
            fundamental.push_back(k);
    }
    fundamental.push_back(N);

    // Relaxed amalgamation. Chains of small supernodes (odometry in pose
    // graphs) cost more in kernel calls than in flops, so a supernode is
    // merged into the next one when it is that one's child in the elimination
    // tree (postorder makes the last child adjacent) and the merged supernode
    // stays small and mostly nonzero. The merged columns take the rows of the
    // parent supernode, a superset of their own (lp(c) \ S ⊆ ext(parent)), so
    // the external rows of a supernode are still lower_pattern(last column).
    // Merge if the supernode gets at most always_cols columns, or at most
    // max_cols columns with at most max_zero_frac explicit zeros.
    constexpr Index always_cols = 6;
    constexpr Index max_cols = 24;
    constexpr double max_zero_frac = 0.5;
    std::vector<Index> size(N), offset(N + 1, 0);
    for (Index k = 0; k < N; ++k)
    {
        size[k] = sp.block_size(order.perm()[k]);
        offset[k + 1] = offset[k] + size[k];
    }
    // Scalar nonzeros of the exact L in the columns of [a, b).
    auto exact_nnz = [&](Index a, Index b)
    {
        double nnz = 0.;
        for (Index k = a; k < b; ++k)
        {
            Index ext = 0;
            for (Index i : lower_pattern_[k])
                ext += size[i];
            nnz += 0.5 * size[k] * (size[k] + 1) + double(size[k]) * ext;
        }
        return nnz;
    };
    supernode_start_.push_back(0);
    double group_nnz = exact_nnz(0, fundamental[1]);
    for (size_t f = 1; f + 1 < fundamental.size(); ++f)
    {
        const Index a = supernode_start_.back(), b = fundamental[f], e = fundamental[f + 1];
        const double nnz = exact_nnz(b, e);
        const double cols = offset[e] - offset[a];
        Index ext = 0;
        for (Index i : lower_pattern_[e - 1])
            ext += size[i];
        const double merged_nnz = 0.5 * cols * (cols + 1) + cols * ext;
        const bool merge = parent_[b - 1] == b &&
                           (cols <= always_cols ||
                            (cols <= max_cols &&
                             merged_nnz - group_nnz - nnz <= max_zero_frac * merged_nnz));
        if (merge)
            group_nnz += nnz;
        else
        {
            supernode_start_.push_back(b);
            group_nnz = nnz;
        }
    }
    supernode_start_.push_back(N); // sentinel
}
