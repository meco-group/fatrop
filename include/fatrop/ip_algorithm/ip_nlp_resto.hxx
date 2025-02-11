//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_nlp_resto_hxx__
#define __fatrop_ip_algorithm_ip_nlp_resto_hxx__
#include "fatrop/common/options.hpp"
#include "fatrop/ip_algorithm/ip_nlp_resto.hpp"
#include "fatrop/nlp/dims.hpp"
#include <cmath>

namespace fatrop
{

    template <typename ProblemType>
    IpNlpResto<ProblemType>::IpNlpResto(const NlpSp &nlp)
        : nlp_(nlp), x_reference_(nlp->nlp_dims().number_of_variables +
                                  nlp->nlp_dims().number_of_ineq_constraints),
          dr_(x_reference_.m())
    {
        const NlpDims &orig_dims = nlp->nlp_dims();
        dims_.number_of_variables = orig_dims.number_of_variables;
        /* this is really the number of slack variables, here we add for the n and p variables */
        dims_.number_of_ineq_constraints =
            orig_dims.number_of_ineq_constraints + 2 * orig_dims.number_of_eq_constraints;
        dims_.number_of_eq_constraints = orig_dims.number_of_eq_constraints;
    }

    template <typename ProblemType> const NlpDims &IpNlpResto<ProblemType>::nlp_dims() const
    {
        return dims_;
    }

    template <typename ProblemType>
    const ProblemDims<ProblemType> &IpNlpResto<ProblemType>::problem_dims() const
    {
        return nlp_->problem_dims();
    }

    template <typename ProblemType>
    Index IpNlpResto<ProblemType>::eval_lag_hess(const ProblemInfo<ProblemType> &info,
                                                 const Scalar objective_scale,
                                                 const VecRealView &primal_x,
                                                 const VecRealView &primal_s,
                                                 const VecRealView &lam, Hessian<ProblemType> &hess)
    {
        return nlp_->eval_lag_hess(info, 0.0, primal_x, primal_s, lam, hess);
    }

    template <typename ProblemType>
    Index IpNlpResto<ProblemType>::eval_constr_jac(const ProblemInfo<ProblemType> &info,
                                                   const VecRealView &primal_x,
                                                   const VecRealView &primal_s,
                                                   Jacobian<ProblemType> &jac)
    {
        return nlp_->eval_constr_jac(info, primal_x, primal_s, jac);
    }

    template <typename ProblemType>
    Index IpNlpResto<ProblemType>::eval_constraint_violation(const ProblemInfo<ProblemType> &info,
                                                             const VecRealView &primal_x,
                                                             const VecRealView &primal_s,
                                                             VecRealView &res)
    {
        Index ret = nlp_->eval_constraint_violation(info, primal_x, primal_s, res);
        // now we add the +p and -n slack variables where applicable
        VecRealView p = primal_s.block(info.number_of_eq_constraints, info.offset_p);
        VecRealView n = primal_s.block(info.number_of_eq_constraints, info.offset_n);
        res = res + if_else(info.constraint_allows_dual_damping, p - n, VecRealScalar(res.m(), 0.));
        return ret;
    }

    template <typename ProblemType>
    Index IpNlpResto<ProblemType>::eval_objective_gradient(const ProblemInfo<ProblemType> &info,
                                                           const Scalar objective_scale,
                                                           const VecRealView &primal_x,
                                                           const VecRealView &primal_s,
                                                           VecRealView &grad_x, VecRealView &grad_s)
    {
        VecRealView dr_x = dr_.block(info.number_of_primal_variables, 0);
        VecRealView x_reference_x = x_reference_.block(info.number_of_primal_variables, 0);
        grad_x = zeta_ * dr_x * (primal_x - x_reference_x);
        VecRealView dr_s = dr_.block(info.number_of_slack_variables, info.offset_slack);
        VecRealView x_reference_s =
            x_reference_.block(info.number_of_slack_variables, info.offset_slack);
        VecRealView grad_ss = grad_s.block(info.number_of_slack_variables, 0);
        VecRealView primal_s_s = primal_s.block(info.number_of_slack_variables, 0);
        grad_ss = zeta_ * dr_s * (primal_s_s - x_reference_s);

        VecRealView gs_p = grad_s.block(info.number_of_eq_constraints, info.offset_p);
        VecRealView gs_n = grad_s.block(info.number_of_eq_constraints, info.offset_n);
        gs_p = if_else(info.constraint_allows_dual_damping, VecRealScalar(gs_p.m(), 1.),
                       VecRealScalar(gs_p.m(), 0.));
        gs_n = if_else(info.constraint_allows_dual_damping, VecRealScalar(gs_n.m(), 1.),
                       VecRealScalar(gs_n.m(), 0.));

        return 0.;
    }

    template <typename ProblemType>
    Index IpNlpResto<ProblemType>::eval_objective(const ProblemInfo<ProblemType> &info,
                                                  const Scalar objective_scale,
                                                  const VecRealView &primal_x,
                                                  const VecRealView &primal_s, Scalar &res)
    {
        res = 0.;
        VecRealView dr_x = dr_.block(info.number_of_primal_variables, 0);
        VecRealView x_reference_x = x_reference_.block(info.number_of_primal_variables, 0);
        res += 0.5 * zeta_ * sumsqr(dr_x * (primal_x - x_reference_x));
        VecRealView dr_s = dr_.block(info.number_of_slack_variables, info.offset_slack);
        VecRealView x_reference_s =
            x_reference_.block(info.number_of_slack_variables, info.offset_slack);
        VecRealView primal_s_s = primal_s.block(info.number_of_slack_variables, 0);
        res += 0.5 * zeta_ * sumsqr(dr_s * (primal_s_s - x_reference_s));
        VecRealView p = primal_s.block(info.number_of_eq_constraints, info.offset_p);
        VecRealView n = primal_s.block(info.number_of_eq_constraints, info.offset_n);
        res += sum(if_else(info.constraint_allows_dual_damping, p + n, VecRealScalar(p.m(), 0.)));
        return 0;
    }

    template <typename ProblemType>
    Index IpNlpResto<ProblemType>::get_bounds(const ProblemInfo<ProblemType> &info,
                                              VecRealView &lower_bounds, VecRealView &upper_bounds)
    {
        // get the bounds for s
        VecRealView lower_bounds_s = lower_bounds.block(info.number_of_slack_variables, 0);
        VecRealView upper_bounds_s = upper_bounds.block(info.number_of_slack_variables, 0);
        nlp_->get_bounds(info, lower_bounds_s, upper_bounds_s);
        // set the lower bounds for the constraints that allow dual damping to zero otherwise -infty
        lower_bounds.block(info.number_of_eq_constraints, info.offset_p) = if_else(
            info.constraint_allows_dual_damping, VecRealScalar(info.number_of_eq_constraints, 0.),
            VecRealScalar(info.number_of_eq_constraints, -std::numeric_limits<Scalar>::infinity()));
        lower_bounds.block(info.number_of_eq_constraints, info.offset_n) = if_else(
            info.constraint_allows_dual_damping, VecRealScalar(info.number_of_eq_constraints, 0.),
            VecRealScalar(info.number_of_eq_constraints, -std::numeric_limits<Scalar>::infinity()));
        // set the upper bounds to infty
        upper_bounds.block(info.number_of_eq_constraints, info.offset_p) =
            std::numeric_limits<Scalar>::infinity();
        upper_bounds.block(info.number_of_eq_constraints, info.offset_n) =
            std::numeric_limits<Scalar>::infinity();

        return 0;
    }

    template <typename ProblemType>
    Index IpNlpResto<ProblemType>::get_initial_primal(const ProblemInfo<ProblemType> &info,
                                                      VecRealView &primal_x)
    {
        return nlp_->get_initial_primal(info, primal_x);
    }

    template <typename ProblemType>
    void IpNlpResto<ProblemType>::get_primal_damping(const ProblemInfo<ProblemType> &info,
                                                     VecRealView &damping)
    {
        VecRealView xs_orig =
            damping.block(info.number_of_primal_variables + info.number_of_slack_variables, 0);
        nlp_->get_primal_damping(info, xs_orig);
        xs_orig = xs_orig + dr_;
        damping.block(info.number_of_eq_constraints, info.offset_slack_p) = 0.;
        damping.block(info.number_of_eq_constraints, info.offset_slack_n) = 0.;
    }

    template <typename ProblemType>
    void IpNlpResto<ProblemType>::set_xs_reference(const ProblemInfo<ProblemType> &info,
                                                   const VecRealView &x_reference,
                                                   const VecRealView &s_reference)
    {
        x_reference_.block(info.number_of_primal_variables, 0) = x_reference;
        x_reference_.block(info.number_of_slack_variables, info.offset_slack) = s_reference;
        dr_ = 1. / max(VecRealScalar(x_reference_.m(), 1.), abs(x_reference_));
    }
    template <typename ProblemType>
    void IpNlpResto<ProblemType>::apply_jacobian_s_transpose(const ProblemInfo<ProblemType> &info,
                                                             const VecRealView &multipliers,
                                                             const Scalar alpha,
                                                             const VecRealView &y, VecRealView &out)
    {
        VecRealView y_orig = y.block(info.number_of_slack_variables, 0);
        VecRealView out_orig = y.block(info.number_of_slack_variables, 0);
        nlp_->apply_jacobian_s_transpose(info, multipliers, alpha, y, out);
        VecRealView y_p = y.block(info.number_of_eq_constraints, info.offset_p);
        VecRealView out_p = out.block(info.number_of_eq_constraints, info.offset_p);
        VecRealView y_n = y.block(info.number_of_eq_constraints, info.offset_n);
        VecRealView out_n = out.block(info.number_of_eq_constraints, info.offset_n);
        out_p = alpha * y_p + if_else(info.constraint_allows_dual_damping, multipliers,
                                      VecRealScalar(y_p.m(), 0.));
        out_n = alpha * y_n - if_else(info.constraint_allows_dual_damping, multipliers,
                                      VecRealScalar(y_n.m(), 0.));
    }

    template <typename ProblemType>
    void IpNlpResto<ProblemType>::register_options(OptionRegistry &registry)
    {
    }

} // namespace fatrop

#endif // __fatrop_ip_algorithm_ip_nlp_resto_hxx__
