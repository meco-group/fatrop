//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_data_hpp__
#define __fatrop_ip_algorithm_ip_data_hpp__
#include "fatrop/context/context.hpp"
#include "fatrop/ip_algorithm/ip_iterate.hpp"
namespace fatrop
{
    template <typename ProblemType> struct IpData
    {
        typedef std::shared_ptr<Nlp<ProblemType>> NlpSp;
        typedef IpIterate<ProblemType> Iterate;
        IpData(const NlpSp &nlp);
        // swich the trial iterate and current iterate into current iterate and reset the trial
        // iterate
        void accept_trial_iterate();
        Iterate &current_iterate() { return *current_iterate_; }
        const Iterate &current_iterate() const { return *current_iterate_; }
        Iterate &trial_iterate() { return *trial_iterate_; }
        const Iterate &trial_iterate() const { return *trial_iterate_; }
        Iterate &stored_iterate() { return *stored_iterate_; }
        const Iterate &stored_iterate() const { return *stored_iterate_; }
        void set_mu(const Scalar mu);
        Scalar mu() const { return mu_; }
        Index iteration_number() const { return iteration_number_; }
        void set_iteration_number(const Index iteration_number)
        {
            iteration_number_ = iteration_number;
        }

    private:
        Scalar mu_;               ///< Barrier value of the NLP.
        Index iteration_number_;  ///< Number of the current iteration.
        Iterate iterate_data_[3]; ///< Data for the three iterates
        Iterate *current_iterate_;
        Iterate *trial_iterate_;
        Iterate *stored_iterate_;
    };

} // namespace fatrop

#endif //__fatrop_ip_algorithm_ip_data_hpp__