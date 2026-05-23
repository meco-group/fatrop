//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_graph_block_sparsity_hpp__
#define __fatrop_graph_block_sparsity_hpp__

#include "fatrop/context/context.hpp"
#include "fatrop/graph/fwd.hpp"

#include <utility>
#include <vector>

namespace fatrop
{
    /**
     * @brief Block-sparse symmetric structure: a vertex-weighted undirected graph.
     *
     * The graph has @c num_blocks vertices, each with a user-prescribed block
     * size @c block_sizes[k]. Each undirected edge @c (i, j) with @c i != j
     * represents an off-diagonal nonzero block; the symmetric counterpart
     * @c (j, i) is implied. All diagonal blocks @c (k, k) are considered
     * structurally present.
     *
     * The adjacency list @c neighbors(k) returns the sorted list of other
     * blocks connected to @c k. The lower-triangular pattern (used to lay out
     * the matrix in column-major form) is exposed via @c column_pattern(j),
     * which lists the row indices @c i >= j with a stored block in column j.
     */
    class BlockSparsityPattern
    {
    public:
        /**
         * @brief Construct from per-block sizes and a list of undirected edges.
         *
         * @param block_sizes Size @c n_k of every block, length @c num_blocks.
         * @param off_diag_edges Pairs @c (i, j) with @c i != j marking
         *                       structurally nonzero off-diagonal blocks.
         *                       Duplicates and unsorted entries are tolerated.
         */
        BlockSparsityPattern(const std::vector<Index> &block_sizes,
                             const std::vector<std::pair<Index, Index>> &off_diag_edges);

        Index num_blocks() const { return num_blocks_; }
        Index total_size() const { return total_size_; }
        Index block_size(Index k) const { return block_sizes_[k]; }
        Index block_offset(Index k) const { return block_offsets_[k]; }

        const std::vector<Index> &block_sizes() const { return block_sizes_; }
        const std::vector<Index> &block_offsets() const { return block_offsets_; }

        /// Sorted list of blocks @c j != k connected to @c k.
        const std::vector<Index> &neighbors(Index k) const { return neighbors_[k]; }

        /// Sorted list of row indices @c i >= j with a structurally stored
        /// lower-triangular block in column @c j (the diagonal index @c j is
        /// always the first entry).
        const std::vector<Index> &column_pattern(Index j) const { return column_pattern_[j]; }

    private:
        Index num_blocks_;
        Index total_size_;
        std::vector<Index> block_sizes_;
        std::vector<Index> block_offsets_;
        std::vector<std::vector<Index>> neighbors_;
        std::vector<std::vector<Index>> column_pattern_;
    };
} // namespace fatrop

#endif // __fatrop_graph_block_sparsity_hpp__
