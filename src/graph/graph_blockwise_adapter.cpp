//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#include "fatrop/graph/graph_blockwise_abstract.hpp"

#include <algorithm>

using namespace fatrop;

GraphBlockwiseAdapter::GraphBlockwiseAdapter(const std::shared_ptr<GraphBlockwiseAbstract> &graph)
    : graph_(graph)
{
    const Index N = graph_->get_num_blocks();
    structure_.block_sizes.resize(N);
    structure_.block_sizes_tan.resize(N);
    structure_.ng_ineq.resize(N);
    offsets_tan_.assign(N + 1, 0);
    offsets_ineq_.assign(N + 1, 0);
    for (Index k = 0; k < N; ++k)
    {
        structure_.block_sizes[k] = graph_->get_block_size(k);
        structure_.block_sizes_tan[k] = graph_->get_block_size_tan(k);
        structure_.ng_ineq[k] = graph_->get_ng_ineq(k);
        offsets_tan_[k + 1] = offsets_tan_[k] + structure_.block_sizes_tan[k];
        offsets_ineq_[k + 1] = offsets_ineq_[k] + structure_.ng_ineq[k];
    }
    auto &edges = structure_.off_diag_edges;
    for (const auto &e : graph_->get_off_diag_edges())
        if (e.first != e.second)
            edges.emplace_back(std::max(e.first, e.second), std::min(e.first, e.second));
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
}

Index GraphBlockwiseAdapter::eval_hessian(const Scalar *objective_scale, const Scalar *x,
                                          const Scalar *mult, GraphHessianView &H)
{
    const Index N = static_cast<Index>(structure_.block_sizes.size());
    for (Index k = 0; k < N; ++k)
    {
        MatRealView blk = H.diag(k);
        Index status = graph_->eval_Hk(k, k, objective_scale, x, mult, &blk.mat());
        if (status != 0)
            return status;
    }
    const auto &edges = structure_.off_diag_edges;
    for (size_t e = 0; e < edges.size(); ++e)
    {
        MatRealView blk = H.off(static_cast<Index>(e));
        Index status = graph_->eval_Hk(edges[e].first, edges[e].second, objective_scale, x,
                                       mult, &blk.mat());
        if (status != 0)
            return status;
    }
    return eval_gradient(objective_scale, x, H.grad());
}

Index GraphBlockwiseAdapter::eval_gradient(const Scalar *objective_scale, const Scalar *x,
                                           Scalar *grad)
{
    const Index N = static_cast<Index>(structure_.block_sizes.size());
    for (Index k = 0; k < N; ++k)
    {
        if (structure_.block_sizes_tan[k] == 0)
            continue;
        Index status = graph_->eval_grad_k(k, objective_scale, x, grad + offsets_tan_[k]);
        if (status != 0)
            return status;
    }
    return 0;
}

Index GraphBlockwiseAdapter::eval_jacobian_ineq(const Scalar *x, GraphJacobianView &J)
{
    const Index N = static_cast<Index>(structure_.block_sizes.size());
    for (Index k = 0; k < N; ++k)
    {
        if (structure_.ng_ineq[k] == 0)
            continue;
        MatRealView blk = J.block(k);
        Index status = graph_->eval_Ggt_ineq_k(k, x, &blk.mat());
        if (status != 0)
            return status;
    }
    return 0;
}

Index GraphBlockwiseAdapter::eval_gineq(const Scalar *x, Scalar *g)
{
    const Index N = static_cast<Index>(structure_.block_sizes.size());
    for (Index k = 0; k < N; ++k)
    {
        if (structure_.ng_ineq[k] == 0)
            continue;
        Index status = graph_->eval_gineq_k(k, x, g + offsets_ineq_[k]);
        if (status != 0)
            return status;
    }
    return 0;
}
