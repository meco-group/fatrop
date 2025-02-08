//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_algorithm_hpp__
#define __fatrop_ip_algorithm_ip_algorithm_hpp__
#include "fatrop/context/context.hpp"
#include "fatrop/ip_algorithm/fwd.hpp"
#include "fatrop/nlp/fwd.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include <memory>

namespace fatrop
{
    /**
     * @brief Enumeration of possible return flags for the interior point solver.
     */
    enum class IpSolverReturnFlag
    {
        Success,                ///< Optimization successful
        MaxIterExceeded,        ///< Maximum number of iterations exceeded
        StopAtTinyStep,         ///< Algorithm stopped due to tiny step size
        StopAtAcceptablePoint,  ///< Algorithm stopped at an acceptable point
        LocalInfeasibility,     ///< Problem is locally infeasible
        FeasiblePointFound,     ///< A feasible point was found
        DivergingIterates,      ///< Iterates are diverging
        ErrorInStepComputation, ///< Error occurred during step computation
        InvalidOption,          ///< An invalid option was provided
        InternalError,          ///< An internal error occurred
        Unknown                 ///< Unknown error or status
    };

    /**
     * @brief Interior Point Algorithm class.
     *
     * This class implements the main interior point algorithm for solving optimization problems.
     *
     * @tparam ProblemType The type of problem being solved.
     */
    template <typename ProblemType> class IpAlgorithm
    {
        typedef std::shared_ptr<IpSearchDirBase> IpSearchDirSp;
        typedef std::shared_ptr<IpLineSearchBase> IpLineSearchSp;
        typedef std::shared_ptr<IpInitializerBase> IpInitializerSp;
        typedef std::shared_ptr<IpMuUpdateBase> IpMuUpdateSp;
        typedef std::shared_ptr<IpEqMultInitializerBase> IpEqMultInitializerSp;
        typedef std::shared_ptr<IpConvergenceCheckBase> IpConvergenceCheckSp;
        typedef std::shared_ptr<IpIterationOutputBase> IpIterationOutputSp;
        typedef std::shared_ptr<IpData<ProblemType>> IpDataSp;

    public:
        /**
         * @brief Construct a new IpAlgorithm object.
         *
         * @param search_dir Search direction calculator
         * @param linesearch Line search algorithm
         * @param initializer Problem initializer
         * @param mu_update Barrier parameter update strategy
         * @param eq_mult_initializer Equality multiplier initializer
         * @param convergence_check Convergence checker
         * @param iteration_output Iteration output handler
         * @param ip_data Interior point algorithm data
         */
        IpAlgorithm(const IpSearchDirSp &search_dir, const IpLineSearchSp &linesearch,
                    const IpInitializerSp &initializer, const IpMuUpdateSp &mu_update,
                    const IpEqMultInitializerSp &eq_mult_initializer,
                    const IpConvergenceCheckSp &convergence_check,
                    const IpIterationOutputSp &iteration_output, const IpDataSp &ip_data);

        /**
         * @brief Reset the algorithm to its initial state.
         */
        void reset();

        /**
         * @brief Run the optimization algorithm.
         *
         * @param is_resto Whether this is a restoration phase (default: false)
         * @return IpSolverReturnFlag The status of the optimization process
         */
        IpSolverReturnFlag optimize(const bool is_resto = false);

        const ProblemInfo<ProblemType> &info() const;

        const VecRealView &solution_primal() const;
        const VecRealView &solution_dual() const;

    private:
        IpSearchDirSp search_dir_;                  ///< Search direction calculator
        IpLineSearchSp linesearch_;                 ///< Line search algorithm
        IpInitializerSp initializer_;               ///< Iterates initializer
        IpMuUpdateSp mu_update_;                    ///< Barrier parameter update strategy
        IpEqMultInitializerSp eq_mult_initializer_; ///< Equality multiplier initializer
        IpConvergenceCheckSp convergence_check_;    ///< Convergence checker
        IpIterationOutputSp iteration_output_;      ///< Iteration output handler
        IpDataSp ip_data_;                          ///< Interior point algorithm data
    };
} // namespace fatrop

#endif //__fatrop_ip_algorithm_ip_algorithm_hpp__
