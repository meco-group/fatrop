//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_data_hxx__
#define __fatrop_ip_algorithm_ip_data_hxx__
#include "fatrop/common/options.hpp"
#include "ip_data.hpp"

namespace fatrop
{

    template <typename ProblemType>
    IpData<ProblemType>::IpData(const NlpSp &nlp)
        : nlp_(nlp), info_(nlp->problem_dims()), iterate_data_{*this, *this, *this},
          current_iterate_(&iterate_data_[0]), trial_iterate_(&iterate_data_[1]),
          stored_iterate_(&iterate_data_[2]),
          hessian_data_{nlp->problem_dims(), nlp->problem_dims()},
          jacobian_data_{nlp->problem_dims(), nlp->problem_dims()},
          lower_bounds_(nlp->nlp_dims().number_of_ineq_constraints),
          upper_bounds_(nlp->nlp_dims().number_of_ineq_constraints),
          lower_bounded_(nlp->nlp_dims().number_of_ineq_constraints),
          upper_bounded_(nlp->nlp_dims().number_of_ineq_constraints),
          single_lower_bounded_(nlp->nlp_dims().number_of_ineq_constraints),
          single_upper_bounded_(nlp->nlp_dims().number_of_ineq_constraints)
    {
        reset();
    }

    template <typename ProblemType> void IpData<ProblemType>::reset()
    {
        // set the bounds
        get_nlp()->get_bounds(info_, lower_bounds_, upper_bounds_);
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
        tiny_step_flag_ = false;
        iteration_number_ = 0;
        // reset associated iterates
        for (auto &iterate : iterate_data_)
        {
            iterate.reset();
        }
        jacobian_curr_ = &jacobian_data_[0];
        hessian_curr_ = &hessian_data_[0];
        jacobian_stored_ = &jacobian_data_[1];
        hessian_stored_ = &hessian_data_[1];
        timing_statistics().reset();
        current_iterate().set_hessian(hessian_curr_);
        current_iterate().set_jacobian(jacobian_curr_);
        trial_iterate().set_hessian(nullptr);
        trial_iterate().set_jacobian(nullptr);
        stored_iterate().set_hessian(nullptr);
        stored_iterate().set_jacobian(nullptr);
        stored_iterate_is_valid_ = false;
    }
    template <typename ProblemType> void IpData<ProblemType>::accept_trial_iterate()
    {
        // switch trial and current iterate
        std::swap(current_iterate_, trial_iterate_);
        trial_iterate_->reset_evaluated_quantities();
        // set the hessian and jacobian pointers
        current_iterate().set_hessian(hessian_curr_);
        current_iterate().set_jacobian(jacobian_curr_);
        trial_iterate().set_hessian(nullptr);
        trial_iterate().set_jacobian(nullptr);
        stored_iterate().set_hessian(nullptr);
        stored_iterate().set_jacobian(nullptr);
    }

    template <typename ProblemType> void IpData<ProblemType>::store_current_iterate()
    {
        // fatrop_assert_msg(!stored_iterate_is_valid_,
        //                   "Trying to store iterate but there is already an iterate stored.");
        // copy current iterate to stored iterate
        stored_iterate() = current_iterate();
        // set the hess and jac pointers of the current iterate to nullptr such that they are not
        // overwritten for certainty
        current_iterate().set_hessian(nullptr);
        current_iterate().set_jacobian(nullptr);
        // swap jacobian and hessian pointers
        std::swap(hessian_curr_, hessian_stored_);
        std::swap(jacobian_curr_, jacobian_stored_);
        stored_iterate_is_valid_ = true;
    }

    template <typename ProblemType> void IpData<ProblemType>::restore_current_iterate()
    {
        fatrop_assert_msg(stored_iterate_is_valid_,
                          "Trying to restore iterate but stored iterate is not valid.");
        // swap current and stored iterate
        std::swap(current_iterate_, stored_iterate_);
        // also swap the hessian and jacobian pointers
        std::swap(hessian_curr_, hessian_stored_);
        std::swap(jacobian_curr_, jacobian_stored_);
        stored_iterate_is_valid_ = false;
    }

    template <typename ProblemType>
    void IpData<ProblemType>::register_options(OptionRegistry &registry)
    {
        registry.register_option("tolerance", &IpData::set_tolerance, this);
    }

} // namespace fatrop

#endif // __fatrop_ip_algorithm_ip_data_hxx__
