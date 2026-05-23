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

    private:
        Index num_blocks_;
        std::vector<Index> parent_;
        std::vector<std::vector<Index>> lower_pattern_;
    };
} // namespace fatrop

#endif // __fatrop_graph_symbolic_factorization_hpp__
