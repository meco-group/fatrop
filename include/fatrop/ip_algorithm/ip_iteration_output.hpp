//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_iteration_output_hpp__
#define __fatrop_ip_iteration_output_hpp__
#include "fatrop/common/options.hpp"
#include "fatrop/context/context.hpp"
#include "fatrop/ip_algorithm/fwd.hpp"
#include "fatrop/ip_algorithm/ip_data.hpp"
#include <memory>

namespace fatrop
{
    /**
     * @brief Base class for iteration output in interior point algorithms.
     */
    class IpIterationOutputBase
    {
    public:
        virtual ~IpIterationOutputBase() = default;

        /**
         * @brief Print the header for the iteration output.
         */
        virtual void print_header() = 0;

        /**
         * @brief Output information about the current iteration.
         */
        virtual void output_current_iteration() = 0;
    };

    /**
     * @brief Concrete implementation of iteration output for a specific problem type.
     *
     * @tparam ProblemType The type of optimization problem being solved.
     */
    template <typename ProblemType> class IpIterationOutput : public IpIterationOutputBase
    {
    public:
        typedef std::shared_ptr<IpData<ProblemType>> IpDataSp;
        typedef IpIterate<ProblemType> IpIterateType;

        /**
         * @brief Construct a new IpIterationOutput object.
         *
         * @param ipdata Shared pointer to the interior point algorithm data.
         */
        IpIterationOutput(const IpDataSp &ipdata);

        void print_header() override;
        void output_current_iteration() override;

    private:
        /**
         * @brief Print detailed information about a single iteration.
         *
         * @param iter Iteration number
         * @param objective Objective function value
         * @param inf_pr Primal infeasibility
         * @param inf_du Dual infeasibility
         * @param lg_mu Log of the barrier parameter
         * @param d_norm Norm of the search direction
         * @param rg Regularization term
         * @param alpha_du Dual step size
         * @param alpha_pr Primal step size
         * @param ls Number of line search iterations
         */
        void print_iteration(Index iter, Scalar objective, Scalar inf_pr, Scalar inf_du,
                             Scalar lg_mu, Scalar d_norm, Scalar rg, Scalar alpha_du,
                             Scalar alpha_pr, Index ls, char alpha_pr_type);

        IpDataSp ipdata_; ///< Shared pointer to the interior point algorithm data

    public:
        // Register options
        void register_options(OptionRegistry &registry);
    };

} // namespace fatrop

#endif // __fatrop_ip_iteration_output_hpp__
