//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_mu_update_hxx__
#define __fatrop_ip_mu_update_hxx__

#include "fatrop/ip_algorithm/ip_data.hpp"
#include "fatrop/ip_algorithm/ip_iterate.hpp"
#include "fatrop/ip_algorithm/ip_linesearch.hpp"
#include "fatrop/ip_algorithm/ip_mu_update.hpp"

namespace fatrop
{
    template <typename ProblemType> void IpMonotoneMuUpdate<ProblemType>::reset()
    {
        ipdata_->current_iterate().set_mu(mu_init_);
        initialized_ = false;
    }

    template <typename ProblemType> bool IpMonotoneMuUpdate<ProblemType>::update_barrier_parameter()
    {
        IpIterateType &curr_it = ipdata_->current_iterate();
        IpIterateType &trial_it = ipdata_->trial_iterate();
        Scalar mu = curr_it.mu();
        Scalar tau = curr_it.tau();
        Scalar sub_problem_error = curr_it.e_mu(mu);
        Scalar kappa_eps_mu = barrier_tol_factor_ * mu;
        Scalar tol = ipdata_->tolerance();
        bool tiny_step_flag = ipdata_->tiny_step_flag();
        ipdata_->set_tiny_step_flag(false);
        bool done = false;
        while ((sub_problem_error < kappa_eps_mu || tiny_step_flag) && !done)
        {
            Scalar new_mu = std::min(mu_linear_decrease_factor_ * mu,
                                     std::pow(mu, mu_superlinear_decrease_power_));
            new_mu = std::max(new_mu, std::min(tol, compl_inf_tol_) / (barrier_tol_factor_ + 1.));
            curr_it.set_mu(new_mu);
            trial_it.set_mu(new_mu);
            bool mu_changed = new_mu != mu;
            mu = new_mu;
            // if first iteration or
            if (initialized_ && !mu_allow_fast_monotone_decrease_)
            {
                done = true;
            }
            else if (!mu_changed)
            {
                done = true;
            }
            else
            {
                sub_problem_error = curr_it.e_mu(mu);
                kappa_eps_mu = barrier_tol_factor_ * mu;
                done = (sub_problem_error > kappa_eps_mu);
            }
            if (done && mu_changed)
            {
                linesearch_->reset_linesearch();
            }
            tiny_step_flag = false;
        }
        trial_it.set_mu(mu);
        initialized_ = true;
        return true;
    }

} // namespace fatrop
#endif // __fatrop_ip_mu_update_hxx__
