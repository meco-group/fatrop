//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_iterate_hxx__
#define __fatrop_ip_algorithm_ip_iterate_hxx__

#include "fatrop/common/options.hpp"
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
        nlp_->get_bounds(info_, lower_bounds_, upper_bounds_);
        // set the bound flags
        number_of_bounds_ = 0;
        for (Index i = 0; i < lower_bounds_.m(); i++)
        {
            lower_bounded_[i] = !std::isinf(lower_bounds_(i));
            upper_bounded_[i] = !std::isinf(upper_bounds_(i));
            if (lower_bounded_[i])
                number_of_bounds_++;
            if (upper_bounded_[i])
                number_of_bounds_++;
            bool single_bounded = lower_bounded_[i] ^ upper_bounded_[i];
            single_lower_bounded_[i] = single_bounded && lower_bounded_[i];
            single_upper_bounded_[i] = single_bounded && upper_bounded_[i];
        }
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
        primal_damping_evaluated_ = false;
    }

    template <typename ProblemType> Scalar IpIterate<ProblemType>::obj_value()
    {
        if (!obj_value_evaluated_)
        {
            Index status =
                nlp_->eval_objective(info_, objective_scale, primal_x_, primal_s_, obj_value_);
            fatrop_assert_msg(status == 0, "Error in evaluating the objective function.");
            obj_value_evaluated_ = true;
        }
        return obj_value_;
    }

    template <typename ProblemType> Scalar IpIterate<ProblemType>::barrier_value()
    {
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
    template <typename ProblemType> const VecRealView &IpIterate<ProblemType>::obj_gradient_x()
    {
        if (!obj_gradient_evaluated_)
        {
            Index status = nlp_->eval_objective_gradient(
                info_, objective_scale, primal_x_, primal_s_, obj_gradient_x_, obj_gradient_s_);
            fatrop_assert_msg(status == 0,
                              "Error in evaluating the gradient of the objective function.");
            obj_gradient_evaluated_ = true;
        }
        return obj_gradient_x_;
    }
    template <typename ProblemType> const VecRealView &IpIterate<ProblemType>::obj_gradient_s()
    {
        if (!obj_gradient_evaluated_)
        {
            Index status = nlp_->eval_objective_gradient(
                info_, objective_scale, primal_x_, primal_s_, obj_gradient_x_, obj_gradient_s_);
            fatrop_assert_msg(status == 0,
                              "Error in evaluating the gradient of the objective function.");
            obj_gradient_evaluated_ = true;
        }
        return obj_gradient_s_;
    }
    template <typename ProblemType> const VecRealView &IpIterate<ProblemType>::constr_viol()
    {
        if (!constr_viol_evaluated_)
        {
            Index status =
                nlp_->eval_constraint_violation(info_, primal_x_, primal_s_, constr_viol_);
            fatrop_assert_msg(status == 0, "Error in evaluating the constraint violation.");
            constr_viol_evaluated_ = true;
        }
        return constr_viol_;
    }
    template <typename ProblemType>
    const VecRealView &IpIterate<ProblemType>::dual_infeasibility_x()
    {
        if (!dual_infeasibility_x_evaluated_)
        {
            dual_infeasibility_x_.block(dual_infeasibility_x_.m(), 0) = obj_gradient_x();
            jacobian().transpose_apply_on_right(info_, dual_eq_, 1.0, dual_infeasibility_x_,
                                                dual_infeasibility_x_);
            dual_infeasibility_x_evaluated_ = true;
        }
        return dual_infeasibility_x_;
    }

    template <typename ProblemType>
    const VecRealView &IpIterate<ProblemType>::dual_infeasibility_s()
    {
        // note that the setters of dual bounds make sure that duals where no bound is present are
        // zero
        if (!dual_infeasibility_s_evaluated_)
        {
            dual_infeasibility_s_ = obj_gradient_s() + dual_bounds_u_ - dual_bounds_l_;
            nlp_->apply_jacobian_s_transpose(info_, dual_eq_, 1.0, dual_infeasibility_s_,
                                             dual_infeasibility_s_);
            dual_infeasibility_s_evaluated_ = true;
        }
        return dual_infeasibility_s_;
    }

    template <typename ProblemType> const VecRealView &IpIterate<ProblemType>::barrier_gradient()
    {
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
    template <typename ProblemType> const VecRealView &IpIterate<ProblemType>::delta_lower()
    {
        if (!delta_lower_evaluated_)
        {
            delta_lower_ = if_else(lower_bounded_, primal_s() - lower_bounds_,
                                   VecRealScalar(lower_bounds_.m(), 1.));
            delta_lower_evaluated_ = true;
        }
        return delta_lower_;
    }
    template <typename ProblemType> const VecRealView &IpIterate<ProblemType>::delta_upper()
    {
        if (!delta_upper_evaluated_)
        {
            delta_upper_ = if_else(upper_bounded_, upper_bounds_ - primal_s(),
                                   VecRealScalar(upper_bounds_.m(), 1.));
            delta_upper_evaluated_ = true;
        }
        return delta_upper_;
    }

    template <typename ProblemType>
    const VecRealView &IpIterate<ProblemType>::relaxed_complementarity_l(const Scalar mu_in)
    {
        relaxed_complementarity_l_ = if_else(lower_bounded(), complementarity_l() - mu_in,
                                             VecRealScalar(dual_bounds_l_.m(), 0.));
        return relaxed_complementarity_l_;
    }

    template <typename ProblemType>
    const VecRealView &IpIterate<ProblemType>::relaxed_complementarity_u(const Scalar mu_in)
    {
        relaxed_complementarity_u_ = if_else(upper_bounded(), complementarity_u() - mu_in,
                                             VecRealScalar(dual_bounds_l_.m(), 0.));
        return relaxed_complementarity_u_;
    }

    template <typename ProblemType>
    const VecRealView &IpIterate<ProblemType>::relaxed_complementarity_l()
    {
        return relaxed_complementarity_l(mu());
    }

    template <typename ProblemType>
    const VecRealView &IpIterate<ProblemType>::relaxed_complementarity_u()
    {
        return relaxed_complementarity_u(mu());
    }

    template <typename ProblemType> const VecRealView &IpIterate<ProblemType>::complementarity_l()
    {
        if (!complementarity_l_evaluated_)
        {
            complementarity_l_ = if_else(lower_bounded(), dual_bounds_l_ * delta_lower(),
                                         VecRealScalar(dual_bounds_l_.m(), 0.));
            complementarity_l_evaluated_ = true;
        }
        return complementarity_l_;
    }

    template <typename ProblemType> const VecRealView &IpIterate<ProblemType>::complementarity_u()
    {
        if (!complementarity_u_evaluated_)
        {
            complementarity_u_ = if_else(upper_bounded(), dual_bounds_u_ * delta_upper(),
                                         VecRealScalar(dual_bounds_u_.m(), 0.));
            complementarity_u_evaluated_ = true;
        }
        return complementarity_u_;
    }

    template <typename ProblemType> void IpIterate<ProblemType>::set_mu(Scalar value)
    {
        barrier_value_evaluated_ = false;
        barrier_gradient_evaluated_ = false;
        linear_decrease_barrier_evaluated_ = false;
        dual_infeasibility_s_evaluated_ = false;
        mu_ = value;
    };

    template <typename ProblemType> const VecRealView &IpIterate<ProblemType>::primal_damping()
    {
        if (!primal_damping_evaluated_)
        {
            nlp_->get_primal_damping(info_, primal_damping_);
            primal_damping_evaluated_ = true;
        }
        return primal_damping_;
    }
    template <typename ProblemType> Scalar IpIterate<ProblemType>::linear_decrease_objective()
    {
        if (!linear_decrease_objective_evaluated_)
        {
            linear_decrease_objective_ =
                dot(obj_gradient_x(), delta_primal_x()) + dot(obj_gradient_s(), delta_primal_s());
            linear_decrease_objective_evaluated_ = true;
        }
        return linear_decrease_objective_;
    }
    template <typename ProblemType> Scalar IpIterate<ProblemType>::linear_decrease_barrier()
    {
        if (!linear_decrease_barrier_evaluated_)
        {
            linear_decrease_barrier_ = dot(barrier_gradient(), delta_primal_s());
            linear_decrease_barrier_evaluated_ = true;
        }
        return linear_decrease_barrier_;
    }

    template <typename ProblemType> Scalar IpIterate<ProblemType>::e_mu(Scalar mu)
    {
        Scalar zl1 = norm_l1(dual_bounds_l()) + norm_l1(dual_bounds_u());
        Index number_of_eq_constraints = nlp()->nlp_dims().number_of_eq_constraints;
        Index number_of_dual_vars = number_of_bounds_ + number_of_eq_constraints;
        Scalar lam_mean = (zl1 + norm_l1(dual_eq())) / (number_of_dual_vars);
        Scalar z_mean = zl1 / number_of_bounds_;
        Scalar constraint_violation = norm_inf(constr_viol());
        Scalar dual_infeasibility_x_linf = norm_inf(dual_infeasibility_x());
        Scalar dual_infeasibility_s_linf = norm_inf(dual_infeasibility_s());
        Scalar dual_infeasibility_linf =
            std::max(dual_infeasibility_x_linf, dual_infeasibility_s_linf);
        Scalar complementarity_l_linf = norm_inf(relaxed_complementarity_l(mu));
        Scalar complementarity_u_linf = norm_inf(relaxed_complementarity_u(mu));
        Scalar complementarity_linf = std::max(complementarity_l_linf, complementarity_u_linf);
        Scalar res = 0.;
        Scalar sd = 0.;
        Scalar sc = 0.;
        if (lam_mean > smax_)
        {
            sd = lam_mean / smax_;
            dual_infeasibility_linf = dual_infeasibility_linf / sd;
        }
        if (z_mean > smax_)
        {
            sc = z_mean / smax_;
            complementarity_linf = complementarity_linf / sc;
        }
        res = std::max({constraint_violation, dual_infeasibility_linf, complementarity_linf});
        return res;
    }

    template <typename ProblemType>
    Scalar IpIterate<ProblemType>::maximum_step_size_primal(const Scalar tau)
    {
        return maximum_step_size_primal(tau, delta_primal_s_);
    }

    template <typename ProblemType>
    Scalar IpIterate<ProblemType>::maximum_step_size_primal(const Scalar tau,
                                                            const VecRealView &delta_s)
    {
        Scalar alpha_max_pr = 1.;
        delta_lower();
        delta_upper();
        for (Index i = 0; i < primal_s_.m(); i++)
        {
            Scalar delta_si = delta_s(i);
            if (lower_bounded_[i])
            {
                if (delta_si < 0)
                    alpha_max_pr = std::min(alpha_max_pr, -tau * delta_lower_(i) / delta_si);
            }
            if (upper_bounded_[i])
            {
                if (delta_si > 0)
                    alpha_max_pr = std::min(alpha_max_pr, tau * delta_upper_(i) / delta_si);
            }
        }
        return alpha_max_pr;
    }

    template <typename ProblemType>
    Scalar IpIterate<ProblemType>::maximum_step_size_dual(const Scalar tau)
    {
        return maximum_step_size_dual(tau, delta_dual_bounds_l_, delta_dual_bounds_u_);
    }

    template <typename ProblemType>
    Scalar IpIterate<ProblemType>::maximum_step_size_dual(const Scalar tau,
                                                          const VecRealView &delta_dual_bounds_l,
                                                          const VecRealView &delta_dual_bounds_u)
    {
        Scalar alpha_max_du = 1.;
        for (Index i = 0; i < primal_s_.m(); i++)
        {
            if (lower_bounded_[i])
            {
                Scalar delta_zi = delta_dual_bounds_l(i);
                if (delta_zi < 0)
                    alpha_max_du = std::min(alpha_max_du, -tau * dual_bounds_l_(i) / delta_zi);
            }
            if (upper_bounded_[i])
            {
                Scalar delta_zi = delta_dual_bounds_u(i);
                if (delta_zi < 0)
                    alpha_max_du = std::min(alpha_max_du, -tau * dual_bounds_u_(i) / delta_zi);
            }
        }
        return alpha_max_du;
    }

    template <typename ProblemType>
    std::pair<Scalar, Scalar> IpIterate<ProblemType>::maximum_step_size(const Scalar tau)
    {
        return {maximum_step_size_primal(tau), maximum_step_size_dual(tau)};
    }
    template <typename ProblemType>
    IpIterate<ProblemType>::IpIterate(const NlpSp &nlp)
        : info_(nlp->problem_dims()), primal_x_(nlp->nlp_dims().number_of_variables),
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
          relaxed_complementarity_l_(nlp->nlp_dims().number_of_ineq_constraints),
          relaxed_complementarity_u_(nlp->nlp_dims().number_of_ineq_constraints),
          primal_damping_(nlp->nlp_dims().number_of_variables +
                          nlp->nlp_dims().number_of_ineq_constraints),
          hessian_(nlp->problem_dims()), lower_bounds_(nlp->nlp_dims().number_of_ineq_constraints),
          upper_bounds_(nlp->nlp_dims().number_of_ineq_constraints),
          lower_bounded_(nlp->nlp_dims().number_of_ineq_constraints),
          upper_bounded_(nlp->nlp_dims().number_of_ineq_constraints),
          single_lower_bounded_(nlp->nlp_dims().number_of_ineq_constraints),
          single_upper_bounded_(nlp->nlp_dims().number_of_ineq_constraints),
          Dx_(nlp->nlp_dims().number_of_variables + nlp->nlp_dims().number_of_ineq_constraints),
          De_(nlp->nlp_dims().number_of_eq_constraints)
    {
        initialize();
    }
    template <typename ProblemType> Hessian<ProblemType> &IpIterate<ProblemType>::hessian()
    {
        if (!hessian_evaluated_)
        {
            Index status = nlp_->eval_lag_hess(info_, objective_scale, primal_x_, primal_s_,
                                               dual_eq_, hessian_);
            fatrop_assert_msg(status == 0, "Error in evaluating the Hessian of the Lagrangian.");
            hessian_evaluated_ = true;
        }
        return hessian_;
    }
    template <typename ProblemType> Hessian<ProblemType> &IpIterate<ProblemType>::zero_hessian()
    {
        primal_damping_ = 0.;
        hessian_.set_zero();
        hessian_evaluated_ = false;
        return hessian_;
    }
    template <typename ProblemType> Jacobian<ProblemType> &IpIterate<ProblemType>::jacobian()
    {
        if (!jacobian_evaluated_)
        {
            Index status = nlp_->eval_constr_jac(info_, primal_x_, primal_s_, jacobian_);
            fatrop_assert_msg(status == 0, "Error in evaluating the Jacobian of the constraints.");
            jacobian_evaluated_ = true;
        }
        return jacobian_;
    }
    template <typename ProblemType>
    void IpIterate<ProblemType>::register_options(OptionRegistry &registry)
    {
        registry.register_option("kappa_d", &IpIterate::set_kappa_d, this);
        registry.register_option("smax", &IpIterate::set_smax, this);
    }

} // namespace fatrop

#endif // __fatrop_algorithm_ip_iterate_hxx__
