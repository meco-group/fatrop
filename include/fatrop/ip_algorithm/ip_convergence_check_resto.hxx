//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_convergence_check_resto_hxx__
#define __fatrop_ip_convergence_check_resto_hxx__

#include "fatrop/common/options.hpp"
#include "fatrop/ip_algorithm/ip_convergence_check_resto.hpp"
#include "fatrop/ip_algorithm/ip_linesearch.hpp"

namespace fatrop
{
    template <typename ProblemType>
    IpConvergenceCheckResto<ProblemType>::IpConvergenceCheckResto(const IpDataSp &data_orig,
                                                                  const IpDataSp &data_resto)
        : Base(data_resto), data_orig_(data_orig), data_resto_(data_resto)
    {
    }
    template <typename ProblemType>
    IpConvergenceStatus IpConvergenceCheckResto<ProblemType>::check_converged()
    {
        IpIterateType &curr_it_resto = data_resto_->current_iterate();
        IpIterateType &curr_it_orig = data_orig_->current_iterate();
        IpIterateType &trial_it_resto = data_resto_->trial_iterate();
        IpIterateType &trial_it_orig = data_orig_->trial_iterate();
        // set the trial point for the original problem
        trial_it_orig.set_primal_x(curr_it_resto.primal_x());
        trial_it_orig.set_primal_s(
            curr_it_resto.primal_s().block(data_orig_->info().number_of_slack_variables, 0));

        successive_resto_iter_++;

        IpConvergenceStatus status;

        Scalar orig_trial_theta = norm_l1(trial_it_orig.constr_viol());
        Scalar orig_curr_theta = norm_l1(curr_it_orig.constr_viol());

        Scalar orig_curr_inf_pr =
            std::max(norm_inf(curr_it_orig.constr_viol()), norm_inf(curr_it_orig.constr_viol()));
        Scalar orig_trial_inf_pr =
            std::max(norm_inf(trial_it_orig.constr_viol()), norm_inf(trial_it_orig.constr_viol()));
        Scalar orig_inf_pr_max = std::max(kappa_resto_ * orig_curr_inf_pr,
                                          std::min(data_orig_->tolerance(), constr_viol_tol_));
        if (kappa_resto_ == 0.)
        {
            orig_inf_pr_max = 0.;
        }
        if (first_resto_iter_)
        {
            status = IpConvergenceStatus::Continue;
        }
        else if (orig_trial_inf_pr > orig_inf_pr_max)
        {
            status = IpConvergenceStatus::Continue;
        }
        else
        {
            Scalar orig_trial_barr = trial_it_orig.obj_value() + trial_it_orig.barrier_value();
            auto line_search_orig = line_search_orig_.lock();
            fatrop_assert_msg(
                line_search_orig,
                "No line search object set for original problem in resto convergence check.");
            if (!line_search_orig->is_acceptable_to_filter(orig_trial_barr, orig_trial_theta))
            {
                status = IpConvergenceStatus::Continue;
            }
            else if (!line_search_orig->is_acceptable_to_current_iterate(orig_trial_barr,
                                                                          orig_trial_theta, true))
            {
                status = IpConvergenceStatus::Continue;
            }
            else
            {
                status = IpConvergenceStatus::Converged;
            }
        }
        if (status == IpConvergenceStatus::Continue)
        {
            status = Base::check_converged();
            if (status != IpConvergenceStatus::Continue)
            {
                status = IpConvergenceStatus::RestoFail;
            }
        }
        first_resto_iter_ = false;

        return status;
    }
    template <typename ProblemType> void IpConvergenceCheckResto<ProblemType>::reset()
    {
        Base::reset();
        first_resto_iter_ = true;
        successive_resto_iter_ = 0;
    }
    template <typename ProblemType>
    void IpConvergenceCheckResto<ProblemType>::register_options(OptionRegistry &registry)
    {
        Base::register_options(registry);
        registry.register_option(
            "required_infeasibility_reduction",
            &IpConvergenceCheckResto<ProblemType>::set_required_infeasibility_reduction, this);
        registry.register_option("constr_viol_tol",
                                 &IpConvergenceCheckResto<ProblemType>::set_constr_viol_tol, this);
    }
} // namespace fatrop

#endif // __fatrop_ip_convergence_check_resto_hxx__