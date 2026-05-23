//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_graph_dims_hpp__
#define __fatrop_graph_dims_hpp__

#include "fatrop/context/context.hpp"
#include "fatrop/graph/block_sparsity.hpp"
#include "fatrop/graph/fwd.hpp"
#include "fatrop/graph/problem_type.hpp"
#include "fatrop/nlp/dims.hpp"

#include <memory>
#include <utility>
#include <vector>

namespace fatrop
{
    /**
     * @brief Dimensions of a graph-structured NLP problem.
     *
     * The primal variable @c x is partitioned into @c num_blocks vector blocks
     * of sizes @c block_sizes[k]. Each block may optionally live on a Lie
     * group / manifold of tangent dimension @c block_sizes_tan[k]; the
     * default is @c block_sizes_tan == block_sizes (Euclidean blocks).
     *
     * The Lagrangian Hessian has a symmetric block-sparse structure described
     * by an undirected adjacency graph (@c off_diag_edges) over the
     * @c num_blocks vertices. The Hessian and the search direction live in
     * tangent space, so the block sparsity pattern is built using the tangent
     * sizes.
     *
     * For each block @c k, @c ng_ineq[k] block-local inequality constraints
     * are attached: each such constraint depends on the variables of block
     * @c k only.
     *
     * @note Graph problems do not support equality constraints by design —
     *       the reduced KKT system stays SPD and is solved by the
     *       block-Cholesky solver tagged with @ref GraphType.
     */
    template <> struct ProblemDims<GraphProblem>
    {
        /// Euclidean variant: tangent block sizes default to primal block sizes.
        ProblemDims(const std::vector<Index> &block_sizes,
                    const std::vector<std::pair<Index, Index>> &off_diag_edges,
                    const std::vector<Index> &ng_ineq);

        /// Lie-group / manifold variant: explicit per-block tangent sizes.
        ProblemDims(const std::vector<Index> &block_sizes,
                    const std::vector<Index> &block_sizes_tan,
                    const std::vector<std::pair<Index, Index>> &off_diag_edges,
                    const std::vector<Index> &ng_ineq);

        void check_problem_dimensions() const;

        const BlockSparsityPattern &sparsity() const { return *sparsity_ptr; }

        std::vector<Index> block_sizes;     ///< Primal-space block sizes.
        std::vector<Index> block_sizes_tan; ///< Tangent-space block sizes.
        std::vector<std::pair<Index, Index>> off_diag_edges;
        std::vector<Index> ng_ineq;

        Index num_blocks = 0;       ///< Number of variable blocks @c N.
        Index nx = 0;               ///< Total primal-variable dimension.
        Index nx_tan = 0;           ///< Total tangent-space (search direction) dimension.
        Index ng = 0;               ///< Number of equality constraints — always 0.
        Index ng_ineq_total = 0;    ///< Total inequality-constraint dimension.

        /// Cumulative primal-variable offsets, length @c num_blocks + 1.
        std::vector<Index> block_offsets;
        /// Cumulative tangent-space offsets, length @c num_blocks + 1.
        std::vector<Index> block_offsets_tan;
        /// Cumulative inequality-constraint offsets, length @c num_blocks + 1.
        std::vector<Index> ng_ineq_offsets;

        /// Block sparsity pattern of the Lagrangian Hessian (in tangent space).
        /// Built once at construction; shared with the Hessian and the linear
        /// solver so they all reference the same canonical instance.
        std::shared_ptr<BlockSparsityPattern> sparsity_ptr;
    };

} // namespace fatrop

#endif // __fatrop_graph_dims_hpp__
