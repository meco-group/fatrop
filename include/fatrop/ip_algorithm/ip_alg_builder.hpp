//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_alg_builder__
#define __fatrop_ip_algorithm_ip_alg_builder__

#include "fatrop/common/fwd.hpp"
#include "fatrop/ip_algorithm/fwd.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/nlp/fwd.hpp"
#include <memory>
#include <optional>

namespace fatrop
{
    /**
     * @brief Builder class for creating and configuring an Interior Point Algorithm.
     *
     * @tparam ProblemType The type of optimization problem being solved.
     */
    template <typename ProblemType> class IpAlgBuilder
    {
    public:
        /**
         * @brief Construct a new IpAlgBuilder object.
         *
         * @param nlp Shared pointer to the Nonlinear Programming problem.
         */
        IpAlgBuilder(const std::shared_ptr<Nlp<ProblemType>> &nlp);

        IpAlgBuilder &with_options_registry(OptionRegistry *options_registry)
        {
            options_registry_ = options_registry;
            return *this;
        }

        /**
         * @brief Create the IpData component.
         * @return Reference to this IpAlgBuilder for method chaining.
         */
        IpAlgBuilder &create_ipdata();

        /**
         * @brief Create the ProblemInfo component.
         * @return Reference to this IpAlgBuilder for method chaining.
         */
        IpAlgBuilder &create_problem_info();

        /**
         * @brief Create the AugSystemSolver component.
         * @return Reference to this IpAlgBuilder for method chaining.
         */
        IpAlgBuilder &create_aug_system_solver();

        /**
         * @brief Create the PdSolverOrig component.
         * @return Reference to this IpAlgBuilder for method chaining.
         */
        IpAlgBuilder &create_pdsolver();

        /**
         * @brief Create the LinearSolver component.
         * @return Reference to this IpAlgBuilder for method chaining.
         */
        IpAlgBuilder &create_linearsolver();

        /**
         * @brief Create the IpSearchDir component.
         * @return Reference to this IpAlgBuilder for method chaining.
         */
        IpAlgBuilder &create_search_dir();

        /**
         * @brief Create the IpLinesearch component.
         * @return Reference to this IpAlgBuilder for method chaining.
         */
        IpAlgBuilder &create_linesearch();

        /**
         * @brief Create the IpInitializer component.
         * @return Reference to this IpAlgBuilder for method chaining.
         */
        IpAlgBuilder &create_initializer();

        /**
         * @brief Create the IpMuUpdate component.
         * @return Reference to this IpAlgBuilder for method chaining.
         */
        IpAlgBuilder &create_mu_update();

        /**
         * @brief Create the IpEqMultInitializer component.
         * @return Reference to this IpAlgBuilder for method chaining.
         */
        IpAlgBuilder &create_eq_mult_initializer();

        /**
         * @brief Create the IpConvergenceCheck component.
         * @return Reference to this IpAlgBuilder for method chaining.
         */
        IpAlgBuilder &create_convergence_check();

        /**
         * @brief Create the IpIterationOutput component.
         * @return Reference to this IpAlgBuilder for method chaining.
         */

        IpAlgBuilder &create_iteration_output();
        /**
         * @brief Build and return the fully constructed IpAlgorithm.
         *
         * This method ensures all necessary components are created before
         * constructing and returning the IpAlgorithm object.
         *
         * @return std::shared_ptr<IpAlgorithm> The fully constructed IpAlgorithm.
         */
        std::shared_ptr<IpAlgorithm<ProblemType>> build();

    private:
        std::shared_ptr<IpNlpOrig<ProblemType>> nlp_orig_;
        std::shared_ptr<IpData<ProblemType>> ipdata_;
        std::shared_ptr<ProblemInfo<ProblemType>> problem_info_;
        std::shared_ptr<AugSystemSolver<ProblemType>> aug_system_solver_;
        std::shared_ptr<PdSolverOrig<ProblemType>> pd_solver_;
        std::shared_ptr<LinearSolver<PdSolverOrig<ProblemType>, PdSystemType<ProblemType>>>
            linear_solver_;
        std::shared_ptr<IpSearchDirBase> search_dir_;
        std::shared_ptr<IpLineSearchBase> linesearch_;
        std::shared_ptr<IpInitializerBase> initializer_;
        std::shared_ptr<IpMuUpdateBase> mu_update_;
        std::shared_ptr<IpEqMultInitializerBase> eq_mult_initializer_;
        std::shared_ptr<IpConvergenceCheckBase> convergence_check_;
        std::shared_ptr<IpIterationOutputBase> iteration_output_;
        OptionRegistry *options_registry_ = nullptr;
    };
} // namespace fatrop

#endif //__fatrop_ip_algorithm_ip_alg_builder__
