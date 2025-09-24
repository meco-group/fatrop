//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_convergence_check_hpp__
#define __fatrop_ip_algorithm_ip_convergence_check_hpp__
#include "fatrop/common/fwd.hpp"
#include "fatrop/ip_algorithm/fwd.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/linear_algebra/vector.hpp"
#include <memory>
namespace fatrop
{
    /**
     * @brief Enumeration of possible convergence statuses for the interior point algorithm.
     */
    enum class IpConvergenceStatus
    {
        Continue,                   ///< Algorithm should continue
        Converged,                  ///< Algorithm has converged to optimal solution
        ConvergedToAcceptablePoint, ///< Algorithm has converged to an acceptable suboptimal point
        MaxIterExceeded,            ///< Maximum number of iterations exceeded
        RestoFail,
        Unassigned                  ///< Convergence status not yet determined
    };

    /**
     * @brief Base class for convergence checking in the interior point algorithm.
     */
    class IpConvergenceCheckBase
    {
    public:
        /**
         * @brief Check if the algorithm has converged.
         * @return The current convergence status.
         */
        virtual IpConvergenceStatus check_converged() = 0;
        virtual void register_options(OptionRegistry &registry) = 0;

        /**
         * @brief Check if the current solution is acceptable.
         * @return True if the solution is acceptable, false otherwise.
         */
        virtual bool check_acceptable() const = 0;

        /**
         * @brief Reset the convergence checker to its initial state.
         */
        virtual void reset() = 0;

    protected:
        virtual ~IpConvergenceCheckBase() = default;
    };

    /**
     * @brief Concrete implementation of convergence checking for a specific problem type.
     * @tparam ProblemType The type of optimization problem being solved.
     */
    template <typename ProblemType> class IpConvergenceCheck : public IpConvergenceCheckBase
    {
        typedef std::shared_ptr<IpData<ProblemType>> IpDataSp;

    public:
        /**
         * @brief Construct a new IpConvergenceCheck object.
         * @param ipdata Shared pointer to the interior point algorithm data.
         */
        IpConvergenceCheck(const IpDataSp &ipdata);

        IpConvergenceStatus check_converged() override;
        void reset() override;
        bool check_acceptable() const override;

    private:
        IpDataSp ipdata_;
        Scalar tol_acceptable_ = 1e-6; ///< Tolerance for acceptable point
        Index acceptable_counter_ = 0; ///< Counter for consecutive acceptable iterations
        Index acceptable_iter_ = 15;   ///< Number of consecutive acceptable iterations required
        Index max_iter_ = 1000;        ///< Maximum number of iterations allowed
        Scalar constr_viol_tol_ = 1e-4; ///< Tolerance for constraint violation

    public:
        // Setter methods for options
        void set_tol_acceptable(const Scalar &value) { tol_acceptable_ = value; }
        void set_acceptable_iter(const Index &value) { acceptable_iter_ = value; }
        void set_max_iter(const Index &value) { max_iter_ = value; }
        void set_constr_viol_tol(const Scalar &value) { constr_viol_tol_ = value; }

        // Register options
        void register_options(OptionRegistry &registry);
    };

} // namespace fatrop

#endif // __fatrop_ip_algorithm_ip_convergence_check_hpp__
