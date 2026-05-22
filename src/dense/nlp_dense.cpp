//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

#include "fatrop/dense/nlp_dense.hpp"
#include "fatrop/dense/hessian.hpp"
#include "fatrop/dense/jacobian.hpp"
#include "fatrop/dense/problem_info.hpp"
#include "fatrop/linear_algebra/linear_algebra.hpp"

using namespace fatrop;

namespace
{
NlpDims make_nlp_dims(const ProblemDims<DenseType> &dense_dims)
{
    return NlpDims(/*number_of_variables=*/dense_dims.nx,
                   /*number_of_tangent_variables=*/dense_dims.nx_tangent,
                   /*number_of_eq_constraints=*/dense_dims.ng + dense_dims.ng_ineq,
                   /*number_of_ineq_constraints=*/dense_dims.ng_ineq);
}
} // namespace

NlpDense::NlpDense(const DenseAbstractSp &dense)
    : dense_(dense),
      dense_dims_(dense->get_nx(), dense->get_nx_tangent(), dense->get_ng(), dense->get_ng_ineq()),
      nlp_dims_(make_nlp_dims(dense_dims_))
{
}

Index NlpDense::eval_lag_hess(const DenseInfo &info, const Scalar objective_scale,
                              const VecRealView &primal_x, const VecRealView & /*primal_s*/,
                              const VecRealView &lam, Hessian<DenseType> &hess)
{
    return dense_->eval_Hh(&objective_scale, primal_x.data(), lam.data(), &hess.Hht.mat());
}

Index NlpDense::eval_constr_jac(const DenseInfo &info, const VecRealView &primal_x,
                                const VecRealView & /*primal_s*/, Jacobian<DenseType> &jac)
{
    dense_->eval_Ggt(primal_x.data(), &jac.Gg_eqt.mat());
    dense_->eval_Ggt_ineq(primal_x.data(), &jac.Gg_ineqt.mat());
    return 0;
}

Index NlpDense::eval_constraint_violation(const DenseInfo &info, const VecRealView &primal_x,
                                          const VecRealView &primal_s, VecRealView &res)
{
    Scalar *res_p = res.data();
    if (info.dims.ng > 0)
        dense_->eval_g(primal_x.data(), res_p);
    if (info.dims.ng_ineq > 0)
        dense_->eval_gineq(primal_x.data(), res_p + info.offset_g_eq_slack);
    // subtract slacks so the slack-equality block is g_ineq(x) - s
    res.block(info.number_of_g_eq_slack, info.offset_g_eq_slack) =
        res.block(info.number_of_g_eq_slack, info.offset_g_eq_slack) -
        primal_s.block(info.number_of_g_eq_slack, 0);
    return 0;
}

Index NlpDense::eval_objective_gradient(const DenseInfo &info, const Scalar objective_scale,
                                        const VecRealView &primal_x,
                                        const VecRealView & /*primal_s*/, VecRealView &grad_x,
                                        VecRealView &grad_s)
{
    grad_s.block(info.number_of_g_eq_slack, 0) = 0;
    return dense_->eval_grad(&objective_scale, primal_x.data(), grad_x.data());
}

Index NlpDense::eval_objective(const DenseInfo &info, const Scalar objective_scale,
                               const VecRealView &primal_x, const VecRealView & /*primal_s*/,
                               Scalar &res)
{
    return dense_->eval_f(&objective_scale, primal_x.data(), &res);
}

Index NlpDense::get_bounds(const DenseInfo &info, VecRealView &lower_bounds,
                           VecRealView &upper_bounds)
{
    if (info.dims.ng_ineq == 0)
        return 0;
    return dense_->get_bounds(lower_bounds.data(), upper_bounds.data());
}

Index NlpDense::get_initial_primal(const DenseInfo &info, VecRealView &primal_x)
{
    return dense_->get_initial(primal_x.data());
}

void NlpDense::get_primal_damping(const DenseInfo &info, VecRealView &damping) { damping = 0; }

void NlpDense::apply_jacobian_s_transpose(const DenseInfo &info, const VecRealView &multipliers,
                                          const Scalar alpha, const VecRealView &y,
                                          VecRealView &out)
{
    // out = alpha * y, then subtract the (slack-equality) multiplier block.
    out = alpha * y;
    out.block(info.number_of_slack_variables, 0) =
        out.block(info.number_of_slack_variables, 0) -
        multipliers.block(info.number_of_slack_variables, info.offset_g_eq_slack);
}

void NlpDense::apply_retraction(const DenseInfo &info, const VecRealView &primal_x,
                                const VecRealView &delta_primal_x, const Scalar alpha,
                                VecRealView &primal_x_next)
{
    dense_->apply_retraction(primal_x.data(), delta_primal_x.data(), alpha, primal_x_next.data());
}

void NlpDense::apply_dual_eq_transformation(const DenseInfo &info, const VecRealView &primal_x,
                                            const VecRealView &dual_eq_in,
                                            VecRealView &dual_eq_out)
{
    dense_->apply_dual_eq_transformation(primal_x.data(), dual_eq_in.data(), dual_eq_out.data());
}
