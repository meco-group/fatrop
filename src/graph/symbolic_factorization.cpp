//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#include "fatrop/graph/symbolic_factorization.hpp"

#include "fatrop/graph/block_sparsity.hpp"
#include "fatrop/graph/elimination_order.hpp"

#include <algorithm>
#include <set>

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
    std::vector<std::set<Index>> work(N);
    for (Index j = 0; j < N; ++j)
    {
        auto &col = work[j];
        for (Index i : A_lower[j])
            col.insert(i);
        for (Index c : children[j])
        {
            // We only inherit rows strictly greater than j; the entry "j"
            // itself is the parent slot and not stored in the pattern.
            for (Index i : work[c])
                if (i > j)
                    col.insert(i);
        }
        lower_pattern_[j].assign(col.begin(), col.end());
        std::sort(lower_pattern_[j].begin(), lower_pattern_[j].end());
    }
}
