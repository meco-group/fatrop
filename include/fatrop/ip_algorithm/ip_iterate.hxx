//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_iterate_hxx__
#define __fatrop_ip_algorithm_ip_iterate_hxx__

#include "fatrop/ip_algorithm/ip_iterate.hpp"
#include <algorithm>
#include <cmath>
namespace fatrop
{
    template <typename ProblemType> void IpIterate<ProblemType>::initialize()
    {
        // reset evaluated quantities
        reset_evaluated_quantities();
        // set the bounds
        nlp_->get_bounds(info, lower_bounds_, upper_bounds_);
        // set the bound flags
        for (Index i = 0; i < lower_bounds_.m(); i++)
        {
            lower_bounded_[i] = isinf(lower_bounds_(i));
            upper_bounded_[i] = isinf(upper_bounds_(i));
            if (lower_bounded_[i])
                number_of_bounds_++;
            if (upper_bounded_[i])
                number_of_bounds_++;
            bool single_bounded = lower_bounded_[i] ^ upper_bounded_[i];
            single_lower_bounded_[i] = single_bounded && lower_bounded_[i];
            single_upper_bounded_[i] = single_bounded && upper_bounded_[i];
        }

        is_initialized_ = true;
    }
    template <typename ProblemType> void IpIterate<ProblemType>::reset_evaluated_quantities()
    {
        obj_value_evaluated_ = false;
        obj_gradient_evaluated_ = false;
        constr_viol_evaluated_ = false;
        dual_infeasibility_x_evaluated_ = false;
        dual_infeasibility_s_evaluated_ = false;
        jacobian_evaluated_ = false;
        hessian_evaluated_ = false;
        barrier_value_evaluated_ = false;
        barrier_gradient_evaluated_ = false;
        delta_lower_evaluated_ = false;
        delta_upper_evaluated_ = false;
        linear_decrease_objective_evaluated_ = false;
        linear_decrease_barrier_evaluated_ = false;
        complementarity_l_evaluated_ = false;
        complementarity_u_evaluated_ = false;
    }

    template <typename ProblemType> Scalar IpIterate<ProblemType>::obj_value()
    {
        fatrop_dbg_assert(is_initialized_);
        if (!obj_value_evaluated_)
        {
            Index status =
                nlp_->eval_objective(info, objective_scale, primal_x_, primal_s_, obj_value_);
            fatrop_assert_msg(status == 0, "Error in evaluating the objective function.");
            obj_value_evaluated_ = true;
        }
        return obj_value_;
    }

    template <typename ProblemType> Scalar IpIterate<ProblemType>::barrier_value()
    {
        fatrop_dbg_assert(is_initialized_);
        if (!barrier_value_evaluated_)
        {
            barrier_value_ = sum(if_else(lower_bounded_, -1. * log(delta_lower()),
                                         VecRealScalar(primal_s_.m(), 0.)) +
                                 if_else(upper_bounded_, -1. * log(delta_upper()),
                                         VecRealScalar(primal_s_.m(), 0.)) +
                                 if_else(single_lower_bounded_, kappa_d_ * delta_lower(),
                                         VecRealScalar(primal_s_.m(), 0.)) +
                                 if_else(single_upper_bounded_, kappa_d_ * delta_upper(),
                                         VecRealScalar(primal_s_.m(), 0.)));
            barrier_value_ *= mu();
        }
        return barrier_value_;
    }
    template <typename ProblemType> VecRealView &IpIterate<ProblemType>::obj_gradient_x()
    {
        fatrop_dbg_assert(is_initialized_);
        if (!obj_gradient_evaluated_)
        {
            Index status = nlp_->eval_objective_gradient(info, objective_scale, primal_x_,
                                                         obj_gradient_x_, obj_gradient_s_);
            fatrop_assert_msg(status == 0,
                              "Error in evaluating the gradient of the objective function.");
            obj_gradient_evaluated_ = true;
        }
        return obj_gradient_x_;
    }
    template <typename ProblemType> VecRealView &IpIterate<ProblemType>::obj_gradient_s()
    {
        fatrop_dbg_assert(is_initialized_);
        if (!obj_gradient_evaluated_)
        {
            Index status = nlp_->eval_objective_gradient(info, objective_scale, primal_x_,
                                                         obj_gradient_x_, obj_gradient_s_);
            fatrop_assert_msg(status == 0,
                              "Error in evaluating the gradient of the objective function.");
            obj_gradient_evaluated_ = true;
        }
        return obj_gradient_s_;
    }
    template <typename ProblemType> VecRealView &IpIterate<ProblemType>::constr_viol()
    {
        fatrop_dbg_assert(is_initialized_);
        if (!constr_viol_evaluated_)
        {
            Index status =
                nlp_->eval_constraint_violation(info, primal_x_, primal_s_, constr_viol_);
            fatrop_assert_msg(status == 0, "Error in evaluating the constraint violation.");
            constr_viol_evaluated_ = true;
        }
        return constr_viol_;
    }
    template <typename ProblemType> VecRealView &IpIterate<ProblemType>::dual_infeasibility_x()
    {
        fatrop_dbg_assert(is_initialized_);
        if (!dual_infeasibility_x_evaluated_)
        {
            dual_infeasibility_x_.block(dual_infeasibility_x_.m(), 0) = obj_gradient_x();
            jacobian().transpose_apply_on_right(info, dual_eq_, 1.0, dual_infeasibility_x_,
                                                dual_infeasibility_x_);
            dual_infeasibility_x_evaluated_ = true;
        }
        return dual_infeasibility_x_;
    }

    template <typename ProblemType> VecRealView &IpIterate<ProblemType>::dual_infeasibility_s()
    {
        fatrop_dbg_assert(is_initialized_);
        if (!dual_infeasibility_s_evaluated_)
        {
            // by convention of the solver dual bounds are zero when not bounded
            dual_infeasibility_s_ =
                -1. * dual_eq_.block(info.number_of_slack_variables, info.offset_g_eq_slack) +
                dual_bounds_u_ - dual_bounds_l_;
            dual_infeasibility_s_evaluated_ = true;
        }
        return dual_infeasibility_s_;
    }

    template <typename ProblemType> VecRealView &IpIterate<ProblemType>::barrier_gradient()
    {
        fatrop_dbg_assert(is_initialized_);
        if (!barrier_gradient_evaluated_)
        {
            const Index m = primal_s_.m();
            barrier_gradient_ =
                if_else(lower_bounded_, -mu() / delta_lower(), VecRealScalar(m, 0.)) +
                if_else(upper_bounded_, +mu() / delta_upper(), VecRealScalar(m, 0.)) +
                if_else(single_lower_bounded_, VecRealScalar(m, mu() * kappa_d_),
                        VecRealScalar(m, 0.)) +
                if_else(single_upper_bounded_, VecRealScalar(m, -mu() * kappa_d_),
                        VecRealScalar(m, 0.));
            barrier_gradient_evaluated_ = true;
        }
        return barrier_gradient_;
    }
    template <typename ProblemType> VecRealView &IpIterate<ProblemType>::delta_lower()
    {
        fatrop_dbg_assert(is_initialized_);
        if (!delta_lower_evaluated_)
        {
            delta_lower_ = if_else(lower_bounded_, delta_primal_s() - lower_bounds_,
                                   VecRealScalar(lower_bounds_.m(), 1.));
            delta_lower_evaluated_ = true;
        }
        return delta_lower_;
    }
    template <typename ProblemType> VecRealView &IpIterate<ProblemType>::delta_upper()
    {
        fatrop_dbg_assert(is_initialized_);
        if (!delta_upper_evaluated_)
        {
            delta_upper_ = if_else(upper_bounded_, upper_bounds_ - delta_primal_s(),
                                   VecRealScalar(upper_bounds_.m(), 1.));
            delta_upper_evaluated_ = true;
        }
        return delta_upper_;
    }

    template <typename ProblemType> VecRealView &IpIterate<ProblemType>::complementarity_l()
    {
        fatrop_dbg_assert(is_initialized_);
        if (!complementarity_l_evaluated_)
        {
            complementarity_l_ = dual_bounds_l_ * delta_lower();
            complementarity_l_evaluated_ = true;
        }
        return complementarity_l_;
    }

    template <typename ProblemType> VecRealView &IpIterate<ProblemType>::complementarity_u()
    {
        fatrop_dbg_assert(is_initialized_);
        if (!complementarity_u_evaluated_)
        {
            complementarity_u_ = dual_bounds_u_ * delta_upper();
            complementarity_u_evaluated_ = true;
        }
        return complementarity_u_;
    }

    template <typename ProblemType> void IpIterate<ProblemType>::set_mu(Scalar value)
    {
        barrier_value_evaluated_ = false;
        barrier_gradient_evaluated_ = false;
        mu_ = value;
    };
    template <typename ProblemType> Scalar IpIterate<ProblemType>::linear_decrease_objective()
    {
        fatrop_dbg_assert(is_initialized_);
        if (!linear_decrease_objective_evaluated_)
        {
            linear_decrease_objective_ = dot(obj_gradient_x(), delta_primal_x());
            +dot(obj_gradient_s(), delta_primal_s());
            linear_decrease_objective_evaluated_ = true;
        }
        return linear_decrease_objective_;
    }
    template <typename ProblemType> Scalar IpIterate<ProblemType>::linear_decrease_barrier()
    {
        fatrop_dbg_assert(is_initialized_);
        if (!linear_decrease_barrier_evaluated_)
        {
            linear_decrease_barrier_ = dot(barrier_gradient(), delta_primal_s());
            linear_decrease_barrier_evaluated_ = true;
        }
        return linear_decrease_barrier_;
    }

    template <typename ProblemType>
    std::pair<Scalar, Scalar> IpIterate<ProblemType>::maximum_step_size(const Scalar tau)
    {
        Scalar alpha_max_pr = 1.;
        Scalar alpha_max_du = 1.;
        for (Index i = 0; i < primal_s_.m(); i++)
        {
            Scalar delta_si = delta_primal_s_(i);
            if (lower_bounded_[i])
            {
                Scalar delta_zi = delta_dual_bounds_l_(i);
                if (delta_si < 0)
                    alpha_max_pr = std::min(alpha_max_pr, -tau * delta_lower_(i) / delta_si);
                if (delta_zi < 0)
                    alpha_max_du = std::min(alpha_max_du, -tau * dual_bounds_l_(i) / delta_zi);
            }
            if (upper_bounded_[i])
            {
                Scalar delta_zi = delta_dual_bounds_u_(i);
                if (delta_si > 0)
                    alpha_max_pr = std::min(alpha_max_pr, tau * delta_upper_(i) / delta_si);
                if (delta_zi < 0)
                    alpha_max_du = std::min(alpha_max_du, -tau * dual_bounds_u_(i) / delta_zi);
            }
        }
        return {alpha_max_pr, alpha_max_du};
    }
    template <typename ProblemType>
    IpIterate<ProblemType>::IpIterate(const NlpSp &nlp)
        : info(nlp->problem_dims()), primal_x_(nlp->nlp_dims().number_of_variables),
          primal_s_(nlp->nlp_dims().number_of_ineq_constraints),
          dual_eq_(nlp->nlp_dims().number_of_eq_constraints),
          dual_bounds_l_(nlp->nlp_dims().number_of_ineq_constraints),
          dual_bounds_u_(nlp->nlp_dims().number_of_ineq_constraints),
          delta_primal_x_(nlp->nlp_dims().number_of_variables),
          delta_primal_s_(nlp->nlp_dims().number_of_ineq_constraints),
          delta_dual_eq_(nlp->nlp_dims().number_of_eq_constraints),
          delta_dual_bounds_l_(nlp->nlp_dims().number_of_ineq_constraints),
          delta_dual_bounds_u_(nlp->nlp_dims().number_of_ineq_constraints),
          obj_gradient_x_(nlp->nlp_dims().number_of_variables),
          obj_gradient_s_(nlp->nlp_dims().number_of_ineq_constraints),
          constr_viol_(nlp->nlp_dims().number_of_eq_constraints),
          dual_infeasibility_x_(nlp->nlp_dims().number_of_variables),
          dual_infeasibility_s_(nlp->nlp_dims().number_of_ineq_constraints),
          barrier_gradient_(nlp->nlp_dims().number_of_ineq_constraints), nlp_(nlp),
          delta_lower_(nlp->nlp_dims().number_of_ineq_constraints),
          delta_upper_(nlp->nlp_dims().number_of_ineq_constraints), jacobian_(nlp->problem_dims()),
          complementarity_l_(nlp->nlp_dims().number_of_ineq_constraints),
          complementarity_u_(nlp->nlp_dims().number_of_ineq_constraints),
          hessian_(nlp->problem_dims()), lower_bounds_(nlp->nlp_dims().number_of_ineq_constraints),
          upper_bounds_(nlp->nlp_dims().number_of_ineq_constraints),
          lower_bounded_(nlp->nlp_dims().number_of_ineq_constraints),
          upper_bounded_(nlp->nlp_dims().number_of_ineq_constraints)
    {
    }
    template <typename ProblemType> Hessian<ProblemType> &IpIterate<ProblemType>::hessian()
    {
        fatrop_dbg_assert(is_initialized_);
        if (!hessian_evaluated_)
        {
            Index status = nlp_->eval_lag_hess(info, objective_scale, primal_x_, primal_s_,
                                               dual_eq_, hessian_);
            fatrop_assert_msg(status == 0, "Error in evaluating the Hessian of the Lagrangian.");
            hessian_evaluated_ = true;
        }
        return hessian_;
    }
    template <typename ProblemType> Jacobian<ProblemType> &IpIterate<ProblemType>::jacobian()
    {
        fatrop_dbg_assert(is_initialized_);
        if (!jacobian_evaluated_)
        {
            Index status = nlp_->eval_constr_jac(info, primal_x_, primal_s_, jacobian_);
            fatrop_assert_msg(status == 0, "Error in evaluating the Jacobian of the constraints.");
            jacobian_evaluated_ = true;
        }
        return jacobian_;
    }
} // namespace fatrop

#endif // __fatrop_algorithm_ip_iterate_hxx__