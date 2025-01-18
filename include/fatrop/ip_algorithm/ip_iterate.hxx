//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_algorithm_ip_iterate_hxx__
#define __fatrop_algorithm_ip_iterate_hxx__

#include "fatrop/ip_algorithm/ip_iterate.hpp"
namespace fatrop
{
    template <typename ProblemType>
    IpIterate<ProblemType>::IpIterate(const NlpSp &nlp)
        : info(nlp->problem_dims()), primal(nlp->nlp_dims().number_of_variables),
          primal_s(nlp->nlp_dims().number_of_ineq_constraints),
          dual_eq(nlp->nlp_dims().number_of_eq_constraints),
          dual_bounds_L(nlp->nlp_dims().number_of_ineq_constraints),
          dual_bounds_Z(nlp->nlp_dims().number_of_ineq_constraints),
          delta_primal(nlp->nlp_dims().number_of_variables),
          delta_primal_s(nlp->nlp_dims().number_of_ineq_constraints),
          delta_dual_eq(nlp->nlp_dims().number_of_eq_constraints),
          delta_dual_bounds_L(nlp->nlp_dims().number_of_ineq_constraints),
          delta_dual_bounds_Z(nlp->nlp_dims().number_of_ineq_constraints),
          obj_gradient(nlp->nlp_dims().number_of_variables),
          constr_viol(nlp->nlp_dims().number_of_variables),
          dual_eq_feas(nlp->nlp_dims().number_of_eq_constraints),
          barrier_gradient(nlp->nlp_dims().number_of_ineq_constraints), nlp_(nlp),
          jacobian_(nlp->problem_dims()), hessian_(nlp->problem_dims())
    {
    }
    template <typename ProblemType> Hessian<ProblemType> &IpIterate<ProblemType>::hessian()
    {
        if (!hessian_evaluated_)
        {
            Index status =
                nlp_->eval_lag_hess(info, objective_scale, primal, primal_s, dual_eq, hessian_);
            fatrop_assert_msg(status == 0, "Error in evaluating the Hessian of the Lagrangian.");
            hessian_evaluated_ = true;
        }
        return hessian_;
    }
    template <typename ProblemType> Jacobian<ProblemType> &IpIterate<ProblemType>::jacobian()
    {
        if (!jacobian_evaluated_)
        {
            Index status = nlp_->eval_constr_jac(info, primal, primal_s, jacobian_);
            fatrop_assert_msg(status == 0, "Error in evaluating the Jacobian of the constraints.");
            jacobian_evaluated_ = true;
        }
        return jacobian_;
    }
    // template <typename ProblemType> VecRealView &IpIterate<ProblemType>::constr_viol()
    // {
    //     if (!constr_viol_evaluated_)
    //     {
    //         Index status = nlp_->eval_constraint_violation(info, primal, primal_s, constr_viol);
    //         fatrop_assert_msg(status == 0, "Error in evaluating the constraint violation.");
    //         constr_viol_evaluated_ = true;
    //     }
    //     return constr_viol;
    // }

} // namespace fatrop

#endif // __fatrop_algorithm_ip_iterate_hxx__