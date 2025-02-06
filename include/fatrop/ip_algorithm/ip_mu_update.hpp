//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_mu_update_hpp__
#define __fatrop_ip_mu_update_hpp__

#include "fatrop/context/context.hpp"
#include "fatrop/ip_algorithm/fwd.hpp"
#include <memory>

namespace fatrop
{
    /**
     * @brief Base class for updating the barrier parameter in interior point methods.
     */
    class IpMuUpdateBase
    {
    public:
        /**
         * @brief Update the barrier parameter.
         *
         * This method should be implemented by derived classes to update
         * the barrier parameter according to the specific update strategy.
         * 
         * @return bool True if the update was successful, false otherwise.
         */
        virtual bool update_barrier_parameter() = 0;

        /**
         * @brief Reset the mu update strategy to its initial state.
         */
        virtual void reset() = 0;

    protected:
        virtual ~IpMuUpdateBase() = default;
    };

    /**
     * @brief Monotone mu (barrier parameter) update strategy for interior point methods.
     *
     * This class implements a monotonically decreasing update strategy for the barrier parameter.
     *
     * @tparam ProblemType The type of optimization problem being solved.
     */
    template <typename ProblemType> class IpMonotoneMuUpdate : public IpMuUpdateBase
    {
        typedef IpData<ProblemType> IpDataType;
        typedef IpIterate<ProblemType> IpIterateType;
        typedef std::shared_ptr<IpDataType> IpDataSp;
        typedef std::shared_ptr<IpLineSearchBase> IpLineSearchSp;

    public:
        /**
         * @brief Construct a new IpMonotoneMuUpdate object.
         *
         * @param ipdata Shared pointer to the interior point algorithm data.
         * @param linesearch Shared pointer to the line search object.
         */
        IpMonotoneMuUpdate(const IpDataSp &ipdata, const IpLineSearchSp &linesearch)
            : ipdata_(ipdata), linesearch_(linesearch)
        {
        }

        /**
         * @brief Update the barrier parameter using a monotone strategy.
         *
         * @return bool True if the update was successful, false otherwise.
         */
        bool update_barrier_parameter() override;

        /**
         * @brief Reset the monotone mu update strategy to its initial state.
         */
        void reset() override;

    private:
        IpDataSp ipdata_;           ///< Shared pointer to the interior point algorithm data.
        IpLineSearchSp linesearch_; ///< Shared pointer to the line search object.
        // optione
        Scalar barrier_tol_factor_ = 10.;
        Scalar mu_linear_decrease_factor_ = 0.2;
        Scalar mu_superlinear_decrease_power_ = 1.5;
        Scalar tau_min_ = 0.99;
        Scalar mu_init_ = 0.1;
        Scalar compl_inf_tol_ = 1e-4;
        bool mu_allow_fast_monotone_decrease_ = true;
        // internal statistics
        bool initialized_ = false;
    };

} // namespace fatrop
#endif // __fatrop_ip_mu_update_hpp__
