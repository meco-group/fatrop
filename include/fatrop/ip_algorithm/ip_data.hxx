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
    }
    template <typename ProblemType> void IpData<ProblemType>::set_mu(const Scalar mu)
    {
        trial_iterate_->set_mu(mu);
        mu_ = mu;
    }
    template <typename ProblemType> void IpData<ProblemType>::accept_trial_iterate()
    {
        Iterate *tmp = &current_iterate();
        current_iterate_ = trial_iterate_;
        trial_iterate_ = tmp;
        trial_iterate_->reset_evaluated_quantities();
    }

} // namespace fatrop

#endif // __fatrop_ip_algorithm_ip_data_hxx__