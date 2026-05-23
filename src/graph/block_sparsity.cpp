//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#include "fatrop/graph/block_sparsity.hpp"

#include "fatrop/common/exception.hpp"

#include <algorithm>

using namespace fatrop;

BlockSparsityPattern::BlockSparsityPattern(
    const std::vector<Index> &block_sizes,
    const std::vector<std::pair<Index, Index>> &off_diag_edges)
    : num_blocks_(static_cast<Index>(block_sizes.size())), total_size_(0),
      block_sizes_(block_sizes), block_offsets_(block_sizes.size() + 1, 0),
      neighbors_(block_sizes.size()), column_pattern_(block_sizes.size())
{
    for (Index k = 0; k < num_blocks_; ++k)
    {
        fatrop_assert(block_sizes_[k] >= 0 && "block sizes must be non-negative");
        block_offsets_[k + 1] = block_offsets_[k] + block_sizes_[k];
    }
    total_size_ = block_offsets_[num_blocks_];

    // Build undirected adjacency (without diagonal self-loops, deduplicated).
    for (const auto &e : off_diag_edges)
    {
        Index a = e.first;
        Index b = e.second;
        fatrop_assert(a >= 0 && a < num_blocks_ && b >= 0 && b < num_blocks_ &&
                      "edge endpoints out of range");
        if (a == b)
            continue;
        neighbors_[a].push_back(b);
        neighbors_[b].push_back(a);
    }
    for (Index k = 0; k < num_blocks_; ++k)
    {
        auto &n = neighbors_[k];
        std::sort(n.begin(), n.end());
        n.erase(std::unique(n.begin(), n.end()), n.end());
    }

    // Build lower-triangular column pattern (diagonal first).
    for (Index j = 0; j < num_blocks_; ++j)
    {
        column_pattern_[j].push_back(j);
        for (Index i : neighbors_[j])
            if (i > j)
                column_pattern_[j].push_back(i);
        std::sort(column_pattern_[j].begin(), column_pattern_[j].end());
    }
}
