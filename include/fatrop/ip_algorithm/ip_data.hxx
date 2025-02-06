//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_data_hxx__
#define __fatrop_ip_algorithm_ip_data_hxx__
#include "ip_data.hpp"

namespace fatrop
{

    template <typename ProblemType>
    IpData<ProblemType>::IpData(const NlpSp &nlp)
        : iterate_data_{nlp, nlp, nlp}, current_iterate_(&iterate_data_[0]),
          trial_iterate_(&iterate_data_[1]), stored_iterate_(&iterate_data_[2])
    {
        reset();
    }

    template <typename ProblemType> void IpData<ProblemType>::reset()
    {
        tiny_step_flag_ = false;
        validate_current_iterate();
    }
    template <typename ProblemType> void IpData<ProblemType>::accept_trial_iterate()
    {
        // switch trial and current iterate
        Iterate *tmp = &current_iterate();
        current_iterate_ = trial_iterate_;
        trial_iterate_ = tmp;
        trial_iterate_->reset_evaluated_quantities();
        validate_current_iterate();
    }

    template <typename ProblemType> void IpData<ProblemType>::backup_current_iterate()
    {
        // switch current and stored iterate
        Iterate *tmp = &current_iterate();
        current_iterate_ = stored_iterate_;
        stored_iterate_ = tmp;
        invalidate_current_iterate();
    }

    template <typename ProblemType> void IpData<ProblemType>::restore_current_iterate()
    {
        // switch current and stored iterate
        Iterate *tmp = &current_iterate();
        current_iterate_ = stored_iterate_;
        stored_iterate_ = tmp;
        validate_current_iterate();
    }

} // namespace fatrop

#endif // __fatrop_ip_algorithm_ip_data_hxx__