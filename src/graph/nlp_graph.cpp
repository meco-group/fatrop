//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#include "fatrop/graph/nlp_graph.hpp"

#include "fatrop/graph/hessian.hpp"
#include "fatrop/graph/jacobian.hpp"
#include "fatrop/graph/problem_info.hpp"
#include "fatrop/linear_algebra/blasfeo_operations.hpp"

using namespace fatrop;

namespace
{
ProblemDims<GraphProblem> make_problem_dims(const GraphAbstract &graph)
{
    const Index N = graph.get_num_blocks();
    std::vector<Index> block_sizes(N);
    std::vector<Index> block_sizes_tan(N);
    std::vector<Index> ng_ineq(N);
    for (Index k = 0; k < N; ++k)
    {
        block_sizes[k] = graph.get_block_size(k);
        block_sizes_tan[k] = graph.get_block_size_tan(k);
        ng_ineq[k] = graph.get_ng_ineq(k);
    }
    return ProblemDims<GraphProblem>(block_sizes, block_sizes_tan, graph.get_off_diag_edges(),
                                     ng_ineq);
}

NlpDims make_nlp_dims(const ProblemDims<GraphProblem> &dims)
{
    return NlpDims(/*number_of_variables=*/dims.nx,
                   /*number_of_tangent_variables=*/dims.nx_tan,
                   /*number_of_eq_constraints=*/dims.ng_ineq_total,
                   /*number_of_ineq_constraints=*/dims.ng_ineq_total);
}
} // namespace

NlpGraph::NlpGraph(const GraphAbstractSp &graph)
    : graph_(graph), dims_(make_problem_dims(*graph)), nlp_dims_(make_nlp_dims(dims_))
{
}

Index NlpGraph::eval_lag_hess(const GraphInfo &info, const Scalar objective_scale,
                              const VecRealView &primal_x, const VecRealView & /*primal_s*/,
                              const VecRealView &lam, Hessian<GraphProblem> &hess)
{
    // Reset Hessian blocks and gradient before writing — the user's eval_Hk is
    // only required to fill its own non-zero block, not zero out neighbours.
    hess.H.set_zero();
    const Scalar *x = primal_x.data();
    const Scalar *mult = lam.data();
    const auto &sp = hess.H.sparsity();
    // Fill structural lower-triangle Hessian blocks.
    for (Index j = 0; j < sp.num_blocks(); ++j)
    {
        for (Index i : sp.column_pattern(j))
        {
            MatRealView blk = hess.H.block(i, j);
            blasfeo_gese_wrap(blk.m(), blk.n(), 0.0, &blk.mat(), blk.ai(), blk.aj());
            Index status = graph_->eval_Hk(i, j, &objective_scale, x, mult, &blk.mat());
            if (status != 0)
                return status;
        }
    }
    // Fill per-block gradient segments (in tangent space).
    for (Index k = 0; k < info.dims.num_blocks; ++k)
    {
        const Index nk_tan = info.dims.block_sizes_tan[k];
        if (nk_tan == 0)
            continue;
        const Index off = info.dims.block_offsets_tan[k];
        Index status = graph_->eval_grad_k(k, &objective_scale, x, hess.grad.data() + off);
        if (status != 0)
            return status;
    }
    return 0;
}

Index NlpGraph::eval_constr_jac(const GraphInfo &info, const VecRealView &primal_x,
                                const VecRealView & /*primal_s*/, Jacobian<GraphProblem> &jac)
{
    const Scalar *x = primal_x.data();
    for (Index k = 0; k < info.dims.num_blocks; ++k)
    {
        if (info.dims.ng_ineq[k] == 0)
            continue;
        MatRealAllocated &blk = jac.Gg_ineqt[k];
        blasfeo_gese_wrap(blk.m(), blk.n(), 0.0, &blk.mat(), 0, 0);
        Index status = graph_->eval_Ggt_ineq_k(k, x, &blk.mat());
        if (status != 0)
            return status;
    }
    return 0;
}

Index NlpGraph::eval_constraint_violation(const GraphInfo &info, const VecRealView &primal_x,
                                          const VecRealView &primal_s, VecRealView &res)
{
    Scalar *res_p = res.data();
    const Scalar *x = primal_x.data();
    for (Index k = 0; k < info.dims.num_blocks; ++k)
    {
        if (info.dims.ng_ineq[k] == 0)
            continue;
        const Index ineq_off = info.offset_g_eq_slack + info.dims.ng_ineq_offsets[k];
        Index status = graph_->eval_gineq_k(k, x, res_p + ineq_off);
        if (status != 0)
            return status;
    }
    // Subtract slacks so the slack-equality block is g_ineq(x) - s.
    res.block(info.number_of_g_eq_slack, info.offset_g_eq_slack) =
        res.block(info.number_of_g_eq_slack, info.offset_g_eq_slack) -
        primal_s.block(info.number_of_g_eq_slack, 0);
    return 0;
}

Index NlpGraph::eval_objective_gradient(const GraphInfo &info, const Scalar objective_scale,
                                        const VecRealView &primal_x,
                                        const VecRealView & /*primal_s*/, VecRealView &grad_x,
                                        VecRealView &grad_s)
{
    grad_s = 0.0;
    const Scalar *x = primal_x.data();
    Scalar *g = grad_x.data();
    for (Index k = 0; k < info.dims.num_blocks; ++k)
    {
        const Index nk_tan = info.dims.block_sizes_tan[k];
        if (nk_tan == 0)
            continue;
        const Index off = info.dims.block_offsets_tan[k];
        Index status = graph_->eval_grad_k(k, &objective_scale, x, g + off);
        if (status != 0)
            return status;
    }
    return 0;
}

Index NlpGraph::eval_objective(const GraphInfo & /*info*/, const Scalar objective_scale,
                                const VecRealView &primal_x, const VecRealView & /*primal_s*/,
                                Scalar &res)
{
    return graph_->eval_f(&objective_scale, primal_x.data(), &res);
}

Index NlpGraph::get_bounds(const GraphInfo &info, VecRealView &lower_bounds,
                            VecRealView &upper_bounds)
{
    if (info.dims.ng_ineq_total == 0)
        return 0;
    return graph_->get_bounds(lower_bounds.data(), upper_bounds.data());
}

Index NlpGraph::get_initial_primal(const GraphInfo & /*info*/, VecRealView &primal_x)
{
    return graph_->get_initial(primal_x.data());
}

void NlpGraph::get_primal_damping(const GraphInfo & /*info*/, VecRealView &damping)
{
    damping = lm_lambda_;
}

void NlpGraph::apply_jacobian_s_transpose(const GraphInfo &info, const VecRealView &multipliers,
                                          const Scalar alpha, const VecRealView &y,
                                          VecRealView &out)
{
    out = alpha * y;
    out.block(info.number_of_slack_variables, 0) =
        out.block(info.number_of_slack_variables, 0) -
        multipliers.block(info.number_of_slack_variables, info.offset_g_eq_slack);
}

void NlpGraph::apply_retraction(const GraphInfo & /*info*/, const VecRealView &primal_x,
                                const VecRealView &delta_primal_x, const Scalar alpha,
                                VecRealView &primal_x_next)
{
    graph_->apply_retraction(primal_x.data(), delta_primal_x.data(), alpha, primal_x_next.data());
}
