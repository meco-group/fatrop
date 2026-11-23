//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_graph_elimination_order_hpp__
#define __fatrop_graph_elimination_order_hpp__

#include "fatrop/context/context.hpp"
#include "fatrop/graph/fwd.hpp"

#include <vector>

namespace fatrop
{
    /**
     * @brief Block-weighted approximate minimum degree elimination order.
     *
     * Computes a permutation of the blocks of @ref BlockSparsityPattern that
     * approximates a minimum-fill ordering for a symmetric block Cholesky
     * factorisation.
     *
     * Runs approximate minimum degree (Amestoy, Davis, Duff 1996) on the
     * quotient graph of the block graph, with degrees weighted by block size:
     * the degree of a block is (an upper bound on) the summed size of its
     * current neighbours, which approximates the cost of eliminating it.
     * Blocks with identical adjacency are merged into supervariables and
     * eliminated consecutively. The result is postordered on the elimination
     * tree, so every subtree occupies consecutive positions.
     *
     * After construction, @c perm()[k] gives the original block index that is
     * eliminated at step @c k, and @c inv_perm()[orig] gives the elimination
     * position of an original block.
     */
    class BlockEliminationOrder
    {
    public:
        explicit BlockEliminationOrder(const BlockSparsityPattern &sp);

        const std::vector<Index> &perm() const { return perm_; }
        const std::vector<Index> &inv_perm() const { return inv_perm_; }

    private:
        std::vector<Index> perm_;
        std::vector<Index> inv_perm_;
    };
} // namespace fatrop

#endif // __fatrop_graph_elimination_order_hpp__
