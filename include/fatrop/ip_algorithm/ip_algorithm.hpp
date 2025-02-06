//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_algorithm_hpp__
#define __fatrop_ip_algorithm_ip_algorithm_hpp__
#include "fatrop/context/context.hpp"
#include "fatrop/ip_algorithm/fwd.hpp"
#include "fatrop/ip_algorithm/ip_iteration_output.hpp"
#include <memory>

namespace fatrop
{
    enum class IpSolverReturnFlag
    {
        Success,
        MaxIterExceeded,
        StopAtTinyStep,
        StopAtAcceptablePoint,
        LocalInfeasibility,
        FeasiblePointFound,
        DivergingIterates,
        ErrorInStepComputation,
        InvalidOption,
        InternalError,
        Unknown
    };
   class IpAlgorithm
    {
        typedef std::shared_ptr<IpSearchDirBase> IpSearchDirSp;
        typedef std::shared_ptr<IpLineSearchBase> IpLineSearchSp;
        typedef std::shared_ptr<IpInitializerBase> IpInitializerSp;
        typedef std::shared_ptr<IpMuUpdateBase> IpMuUpdateSp;
        typedef std::shared_ptr<IpEqMultInitializerBase> IpEqMultInitializerSp;
        typedef std::shared_ptr<IpConvergenceCheckBase> IpConvergenceCheckSp;
        typedef std::shared_ptr<IpIterationOutputBase> IpIterationOutputSp;

    public:
        IpAlgorithm(const IpSearchDirSp &search_dir, const IpLineSearchSp &linesearch,
                    const IpInitializerSp &initializer, const IpMuUpdateSp &mu_update,
                    const IpEqMultInitializerSp &eq_mult_initializer,
                    const IpConvergenceCheckSp &convergence_check,
                    const IpIterationOutputSp &iteration_output);

        void reset();
        IpSolverReturnFlag optimize(const bool is_resto = false);

    private:
        IpSearchDirSp search_dir_;
        IpLineSearchSp linesearch_;
        IpInitializerSp initializer_;
        IpMuUpdateSp mu_update_;
        IpEqMultInitializerSp eq_mult_initializer_;
        IpConvergenceCheckSp convergence_check_;
        IpIterationOutputSp iteration_output_;
    };
} // namespace fatrop

#endif //__fatrop_ip_algorithm_ip_algorithm_hpp__
