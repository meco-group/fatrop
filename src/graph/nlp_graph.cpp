//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#include "fatrop/graph/nlp_graph.hpp"

#include "fatrop/graph/hessian.hpp"
#include "fatrop/graph/jacobian.hpp"
#include "fatrop/graph/problem_info.hpp"
#include "fatrop/common/exception.hpp"

using namespace fatrop;

namespace
{
ProblemDims<GraphProblem> make_problem_dims(const GraphStructure &st)
{
    const Index N = static_cast<Index>(st.block_sizes.size());
    const std::vector<Index> &block_sizes_tan =
        st.block_sizes_tan.empty() ? st.block_sizes : st.block_sizes_tan;
    const std::vector<Index> ng_ineq = st.ng_ineq.empty() ? std::vector<Index>(N, 0) : st.ng_ineq;
    return ProblemDims<GraphProblem>(st.block_sizes, block_sizes_tan, st.off_diag_edges, ng_ineq);
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
    : graph_(graph), dims_(make_problem_dims(graph->get_structure())),
      nlp_dims_(make_nlp_dims(dims_)), diag_ids_(dims_.num_blocks),
      off_ids_(dims_.off_diag_edges.size())
{
    const BlockSparsityPattern &sp = dims_.sparsity();
    for (Index k = 0; k < dims_.num_blocks; ++k)
        diag_ids_[k] = sp.block_id(k, k);
    for (size_t e = 0; e < off_ids_.size(); ++e)
    {
        const auto &edge = dims_.off_diag_edges[e];
        fatrop_assert_msg(edge.first > edge.second, "off_diag_edges must have i > j");
        off_ids_[e] = sp.block_id(edge.first, edge.second);
    }
    fatrop_assert_msg(sp.num_stored_blocks() ==
                          dims_.num_blocks + static_cast<Index>(off_ids_.size()),
                      "off_diag_edges must be unique");
}

Index NlpGraph::eval_lag_hess(const GraphInfo & /*info*/, const Scalar objective_scale,
                              const VecRealView &primal_x, const VecRealView & /*primal_s*/,
                              const VecRealView &lam, Hessian<GraphProblem> &hess)
{
    hess.set_zero();
    GraphHessianView H(hess.H, diag_ids_.data(), off_ids_.data(), hess.grad.data());
    return graph_->eval_hessian(&objective_scale, primal_x.data(), lam.data(), H);
}

Index NlpGraph::eval_constr_jac(const GraphInfo &info, const VecRealView &primal_x,
                                const VecRealView & /*primal_s*/, Jacobian<GraphProblem> &jac)
{
    if (info.dims.ng_ineq_total == 0)
        return 0;
    for (auto &blk : jac.Gg_ineqt)
        blk = 0.0;
    GraphJacobianView J(jac.Gg_ineqt.data());
    return graph_->eval_jacobian_ineq(primal_x.data(), J);
}

Index NlpGraph::eval_constraint_violation(const GraphInfo &info, const VecRealView &primal_x,
                                          const VecRealView &primal_s, VecRealView &res)
{
    if (info.dims.ng_ineq_total > 0)
    {
        Index status = graph_->eval_gineq(primal_x.data(), res.data() + info.offset_g_eq_slack);
        if (status != 0)
            return status;
    }
    // Subtract slacks so the slack-equality block is g_ineq(x) - s.
    res.block(info.number_of_g_eq_slack, info.offset_g_eq_slack) =
        res.block(info.number_of_g_eq_slack, info.offset_g_eq_slack) -
        primal_s.block(info.number_of_g_eq_slack, 0);
    return 0;
}

Index NlpGraph::eval_objective_gradient(const GraphInfo & /*info*/, const Scalar objective_scale,
                                        const VecRealView &primal_x,
                                        const VecRealView & /*primal_s*/, VecRealView &grad_x,
                                        VecRealView &grad_s)
{
    grad_s = 0.0;
    grad_x = 0.0;
    return graph_->eval_gradient(&objective_scale, primal_x.data(), grad_x.data());
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
