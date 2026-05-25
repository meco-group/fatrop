//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_graph_symbolic_factorization_hpp__
#define __fatrop_graph_symbolic_factorization_hpp__

#include "fatrop/context/context.hpp"
#include "fatrop/graph/fwd.hpp"

#include <vector>

namespace fatrop
{
    /**
     * @brief Symbolic factorisation of a block-sparse SPD matrix under a given
     *        elimination order.
     *
     * Given the unpermuted @ref BlockSparsityPattern and an
     * @ref BlockEliminationOrder, computes:
     *   - the elimination tree @c parent[k] (with @c -1 for tree roots), where
     *     indices refer to positions in the elimination order;
     *   - for each column @c j in the elimination order, the sorted list of
     *     row indices @c i > j with a (possibly filled-in) structural nonzero
     *     in the Cholesky factor @c L.
     *
     * All indices in this object are in *permuted* (elimination-order)
     * coordinates. The diagonal index @c j itself is NOT included in
     * @c lower_pattern(j); it is implicit.
     */
    class BlockSymbolicFactorization
    {
    public:
        BlockSymbolicFactorization(const BlockSparsityPattern &sp,
                                   const BlockEliminationOrder &order);

        Index num_blocks() const { return num_blocks_; }

        const std::vector<Index> &parent() const { return parent_; }

        /// Sorted list of row indices @c i > j with a nonzero (incl. fill-in)
        /// in column @c j of @c L, in permuted coordinates.
        const std::vector<Index> &lower_pattern(Index j) const { return lower_pattern_[j]; }

        /// Fundamental supernode partition of the columns.
        ///
        /// A range @c [supernode_start(s), supernode_start(s+1)) of consecutive
        /// columns (in elimination order) forms a supernode iff:
        ///   * each column's elimination-tree parent is the next column, AND
        ///   * the lower pattern of column @c k+1 equals the lower pattern of
        ///     column @c k with the leading row @c k+1 removed
        ///     (equivalently, |lp(k+1)| = |lp(k)| - 1 — the count condition
        ///     suffices once the parent condition holds, because inheritance
        ///     forces @c lp(k+1) ⊇ lp(k) \ {k+1}).
        ///
        /// Columns inside the same supernode share an identical block of
        /// external off-diagonal rows: eliminating them can be done as one
        /// dense (block_size_total + ext) x block_size_total panel factorisation
        /// instead of one panel per column, amortising BLASFEO call overhead.
        /// The sentinel @c supernode_start(num_supernodes()) equals @c num_blocks().
        Index num_supernodes() const { return static_cast<Index>(supernode_start_.size()) - 1; }
        Index supernode_start(Index s) const { return supernode_start_[s]; }
        const std::vector<Index> &supernode_start() const { return supernode_start_; }

    private:
        Index num_blocks_;
        std::vector<Index> parent_;
        std::vector<std::vector<Index>> lower_pattern_;
        std::vector<Index> supernode_start_;
    };
} // namespace fatrop

#endif // __fatrop_graph_symbolic_factorization_hpp__
