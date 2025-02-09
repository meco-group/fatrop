
//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_convergence_check_hxx__
#define __fatrop_ip_algorithm_ip_convergence_check_hxx__

#include "fatrop/ip_algorithm/ip_convergence_check.hpp"
#include "fatrop/ip_algorithm/ip_data.hpp"
#include "fatrop/linear_algebra/vector.hpp"
#include "fatrop/nlp/dims.hpp"
#include "fatrop/ocp/nlp_ocp.hpp"
#include "fatrop/common/options.hpp"

namespace fatrop
{
    template <typename ProblemType>
    IpConvergenceCheck<ProblemType>::IpConvergenceCheck(const IpDataSp &ipdata) : ipdata_(ipdata)
    {
    }

    template <typename ProblemType> void IpConvergenceCheck<ProblemType>::reset()
    {
        acceptable_counter_ = 0;
    }

    template <typename ProblemType> bool IpConvergenceCheck<ProblemType>::check_acceptable() const
    {
        return (ipdata_->current_iterate().e_mu(0.) <= tol_acceptable_);
    }

    template <typename ProblemType>
    IpConvergenceStatus IpConvergenceCheck<ProblemType>::check_converged()
    {
        Scalar tol = ipdata_->tolerance();

        if (ipdata_->current_iterate().e_mu(0.) <= tol)
            return IpConvergenceStatus::Converged;
        if (check_acceptable())
        {
            acceptable_counter_++;
            if (acceptable_counter_ >= acceptable_iter_)
                return IpConvergenceStatus::ConvergedToAcceptablePoint;
        }
        else
        {
            acceptable_counter_ = 0;
        }
        if (ipdata_->iteration_number() >= max_iter_)
            return IpConvergenceStatus::MaxIterExceeded;
        return IpConvergenceStatus::Continue;
    }

    template <typename ProblemType>
    void IpConvergenceCheck<ProblemType>::register_options(OptionRegistry& registry)
    {
        registry.register_option("tol_acceptable", &IpConvergenceCheck::set_tol_acceptable, this);
        registry.register_option("acceptable_iter", &IpConvergenceCheck::set_acceptable_iter, this);
        registry.register_option("max_iter", &IpConvergenceCheck::set_max_iter, this);
    }

} // namespace fatrop
#endif // __fatrop_ip_algorithm_ip_convergence_check_hxx__
