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
     * @brief Block-weighted minimum-degree elimination order.
     *
     * Computes a permutation of the blocks of @ref BlockSparsityPattern that
     * approximates a minimum-fill ordering for a symmetric block Cholesky
     * factorisation.
     *
     * The algorithm runs the classical minimum-degree heuristic on the
     * (quotient) elimination graph, but weights vertices by their block size:
     * the "degree" of a vertex is the sum of the sizes of its current
     * neighbours, which approximates the dense floating-point cost of
     * eliminating that vertex. Eliminating a vertex turns its neighbour set
     * into a clique (modelling block fill-in), then the vertex is removed.
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
