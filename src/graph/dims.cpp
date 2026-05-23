//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#include "fatrop/graph/dims.hpp"

#include "fatrop/common/exception.hpp"

using namespace fatrop;

namespace
{
std::vector<Index> prefix_sum(const std::vector<Index> &sizes)
{
    std::vector<Index> off(sizes.size() + 1, 0);
    for (size_t k = 0; k < sizes.size(); ++k)
        off[k + 1] = off[k] + sizes[k];
    return off;
}
} // namespace

ProblemDims<GraphProblem>::ProblemDims(
    const std::vector<Index> &block_sizes_in,
    const std::vector<std::pair<Index, Index>> &off_diag_edges_in,
    const std::vector<Index> &ng_ineq_in)
    : ProblemDims(block_sizes_in, block_sizes_in, off_diag_edges_in, ng_ineq_in)
{
}

ProblemDims<GraphProblem>::ProblemDims(
    const std::vector<Index> &block_sizes_in, const std::vector<Index> &block_sizes_tan_in,
    const std::vector<std::pair<Index, Index>> &off_diag_edges_in,
    const std::vector<Index> &ng_ineq_in)
    : block_sizes(block_sizes_in), block_sizes_tan(block_sizes_tan_in),
      off_diag_edges(off_diag_edges_in), ng_ineq(ng_ineq_in)
{
    num_blocks = static_cast<Index>(block_sizes.size());
    fatrop_assert_msg(static_cast<Index>(block_sizes_tan.size()) == num_blocks,
                      "block_sizes_tan must have one entry per block");
    fatrop_assert_msg(static_cast<Index>(ng_ineq.size()) == num_blocks,
                      "ng_ineq must have one entry per block");
    block_offsets = prefix_sum(block_sizes);
    block_offsets_tan = prefix_sum(block_sizes_tan);
    ng_ineq_offsets = prefix_sum(ng_ineq);
    nx = block_offsets.back();
    nx_tan = block_offsets_tan.back();
    ng_ineq_total = ng_ineq_offsets.back();
    ng = 0;
    // The sparsity pattern lives in tangent space, so its vertex sizes are
    // the tangent block sizes.
    sparsity_ptr = std::make_shared<BlockSparsityPattern>(block_sizes_tan, off_diag_edges);
    check_problem_dimensions();
}

void ProblemDims<GraphProblem>::check_problem_dimensions() const
{
    for (Index k = 0; k < num_blocks; ++k)
    {
        fatrop_assert_msg(block_sizes[k] >= 0, "block_sizes must be non-negative");
        fatrop_assert_msg(block_sizes_tan[k] >= 0, "block_sizes_tan must be non-negative");
        fatrop_assert_msg(ng_ineq[k] >= 0, "ng_ineq must be non-negative");
    }
    for (const auto &e : off_diag_edges)
    {
        fatrop_assert_msg(e.first != e.second, "off_diag_edges must have distinct endpoints");
        fatrop_assert_msg(e.first >= 0 && e.first < num_blocks, "edge endpoint out of range");
        fatrop_assert_msg(e.second >= 0 && e.second < num_blocks, "edge endpoint out of range");
    }
}
