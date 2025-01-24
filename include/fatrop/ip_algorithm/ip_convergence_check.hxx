
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

namespace fatrop
{
    template <typename ProblemType>
    IpConvergenceCheck<ProblemType>::IpConvergenceCheck(const IpDataSp &ipdata) : ipdata_(ipdata)
    {
    }

    template <typename ProblemType> bool IpConvergenceCheck<ProblemType>::converged() const
    {
        return (ipdata_->current_iterate().e_mu(0.) <= tol_);
    }

} // namespace fatrop
#endif // __fatrop_ip_algorithm_ip_convergence_check_hxx__