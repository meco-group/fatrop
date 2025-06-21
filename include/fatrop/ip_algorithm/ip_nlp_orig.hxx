//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_nlp_orig_hxx__
#define __fatrop_ip_algorithm_ip_nlp_orig_hxx__
#include "fatrop/common/options.hpp"
#include "fatrop/ip_algorithm/ip_nlp_orig.hpp"
#include "fatrop/ip_algorithm/ip_timings.hpp"
#include "fatrop/nlp/dims.hpp"
#include <cmath>

namespace fatrop
{

    template <typename ProblemType>
    IpNlpOrig<ProblemType>::IpNlpOrig(const NlpSp &nlp)
        : nlp_(nlp), modified_bounds_lower_(nlp_->nlp_dims().number_of_ineq_constraints),
          modified_bounds_upper_(nlp_->nlp_dims().number_of_ineq_constraints)
    {
    }

    template <typename ProblemType> const NlpDims &IpNlpOrig<ProblemType>::nlp_dims() const
    {
        return nlp_->nlp_dims();
    }

    template <typename ProblemType>
    const ProblemDims<ProblemType> &IpNlpOrig<ProblemType>::problem_dims() const
    {
        return nlp_->problem_dims();
    }

    template <typename ProblemType>
    Index IpNlpOrig<ProblemType>::eval_lag_hess(const ProblemInfo<ProblemType> &info,
                                                const Scalar objective_scale,
                                                const VecRealView &primal_x,
                                                const VecRealView &primal_s, const VecRealView &lam,
                                                Hessian<ProblemType> &hess)
    {
        if (timings_)
            timings_->eval_hessian.start();
        Index ret = nlp_->eval_lag_hess(info, objective_scale, primal_x, primal_s, lam, hess);
        if (timings_)
            timings_->eval_hessian.pause();
        return ret;
    }

    template <typename ProblemType>
    Index IpNlpOrig<ProblemType>::eval_constr_jac(const ProblemInfo<ProblemType> &info,
                                                  const VecRealView &primal_x,
                                                  const VecRealView &primal_s,
                                                  Jacobian<ProblemType> &jac)
    {
        if (timings_)
            timings_->eval_jacobian.start();
        Index ret = nlp_->eval_constr_jac(info, primal_x, primal_s, jac);
        if (timings_)
            timings_->eval_jacobian.pause();
        return ret;
    }

    template <typename ProblemType>
    Index IpNlpOrig<ProblemType>::eval_constraint_violation(const ProblemInfo<ProblemType> &info,
                                                            const VecRealView &primal_x,
                                                            const VecRealView &primal_s,
                                                            VecRealView &res)
    {
        if (timings_)
            timings_->eval_constraint_violation.start();
        Index ret = nlp_->eval_constraint_violation(info, primal_x, primal_s, res);
        if (timings_)
            timings_->eval_constraint_violation.pause();
        return ret;
    }

    template <typename ProblemType>
    Index IpNlpOrig<ProblemType>::eval_objective_gradient(const ProblemInfo<ProblemType> &info,
                                                          const Scalar objective_scale,
                                                          const VecRealView &primal_x,
                                                          const VecRealView &primal_s,
                                                          VecRealView &grad_x, VecRealView &grad_s)
    {
        if (timings_)
            timings_->eval_gradient.start();
        Index ret = nlp_->eval_objective_gradient(info, objective_scale, primal_x, primal_s, grad_x,
                                                  grad_s);
        if (timings_)
            timings_->eval_gradient.pause();
        return ret;
    }

    template <typename ProblemType>
    Index IpNlpOrig<ProblemType>::eval_objective(const ProblemInfo<ProblemType> &info,
                                                 const Scalar objective_scale,
                                                 const VecRealView &primal_x,
                                                 const VecRealView &primal_s, Scalar &res)
    {
        if (timings_)
            timings_->eval_objective.start();
        Index ret = nlp_->eval_objective(info, objective_scale, primal_x, primal_s, res);
        if (timings_)
            timings_->eval_objective.pause();
        return ret;
    }

    template <typename ProblemType>
    Index IpNlpOrig<ProblemType>::get_bounds(const ProblemInfo<ProblemType> &info,
                                             VecRealView &lower_bounds, VecRealView &upper_bounds)
    {
        const Index number_of_slacks = nlp_->nlp_dims().number_of_ineq_constraints;
        nlp_->get_bounds(info, modified_bounds_lower_, modified_bounds_upper_);
        auto lower_bounded = [&](const Index i) { return !std::isinf(modified_bounds_lower_(i)); };
        auto upper_bounded = [&](const Index i) { return !std::isinf(modified_bounds_upper_(i)); };
        auto updated_lower = modified_bounds_lower_ -
                             min(VecRealScalar(number_of_slacks, constr_viol_tol_),
                                 bound_relax_factor_ * max(VecRealScalar(number_of_slacks, 1.),
                                                           abs(modified_bounds_lower_)));
        auto updated_upper = modified_bounds_upper_ +
                             min(VecRealScalar(number_of_slacks, constr_viol_tol_),
                                 bound_relax_factor_ * max(VecRealScalar(number_of_slacks, 1.),
                                                           abs(modified_bounds_upper_)));
        lower_bounds = updated_lower;
        upper_bounds = updated_upper;
        return 0;
    }

    template <typename ProblemType>
    Index IpNlpOrig<ProblemType>::get_initial_primal(const ProblemInfo<ProblemType> &info,
                                                     VecRealView &primal_x)
    {
        return nlp_->get_initial_primal(info, primal_x);
    }

    template <typename ProblemType>
    void IpNlpOrig<ProblemType>::get_primal_damping(const ProblemInfo<ProblemType> &info,
                                                    VecRealView &damping)
    {
        nlp_->get_primal_damping(info, damping);
    }
    template <typename ProblemType>
    void IpNlpOrig<ProblemType>::apply_jacobian_s_transpose(const ProblemInfo<ProblemType> &info,
                                                            const VecRealView &multipliers,
                                                            const Scalar alpha,
                                                            const VecRealView &y, VecRealView &out)
    {
        nlp_->apply_jacobian_s_transpose(info, multipliers, alpha, y, out);
    }

    template <typename ProblemType>
    void IpNlpOrig<ProblemType>::register_options(OptionRegistry &registry)
    {
        registry.register_option("constr_viol_tol", &IpNlpOrig::set_constr_viol_tol, this);
        registry.register_option("bound_relax_factor", &IpNlpOrig::set_bound_relax_factor, this);
    }

} // namespace fatrop

#endif // __fatrop_ip_algorithm_ip_nlp_orig_hxx__
