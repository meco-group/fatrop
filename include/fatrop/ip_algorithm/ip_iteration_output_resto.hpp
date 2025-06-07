//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_iteration_output_resto_hpp__
#define __fatrop_ip_iteration_output_resto_hpp__
#include "fatrop/common/options.hpp"
#include "fatrop/context/context.hpp"
#include "fatrop/ip_algorithm/fwd.hpp"
#include "fatrop/ip_algorithm/ip_data.hpp"
#include "fatrop/ip_algorithm/ip_iteration_output.hpp"
#include <memory>

namespace fatrop
{
    /**
     * @brief Concrete implementation of iteration output for a specific problem type.
     *
     * @tparam ProblemType The type of optimization problem being solved.
     */
    template <typename ProblemType> class IpIterationOutputResto : public IpIterationOutputBase
    {
    public:
        typedef std::shared_ptr<IpData<ProblemType>> IpDataSp;
        typedef IpIterate<ProblemType> IpIterateType;

        /**
         * @brief Construct a new IpIterationOutputResto object.
         *
         * @param ipdata Shared pointer to the interior point algorithm data.
         */
        IpIterationOutputResto(const IpDataSp &ipdata_orig, const IpDataSp &ipdata_resto);

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

        IpDataSp ipdata_orig_;  ///< Shared pointer to the interior point algorithm data
        IpDataSp ipdata_resto_; ///< Shared pointer to the interior point algorithm data

    public:
        // Register options
        void register_options(OptionRegistry &registry);
    };

} // namespace fatrop

#endif // __fatrop_ip_iteration_output_resto_hpp__
