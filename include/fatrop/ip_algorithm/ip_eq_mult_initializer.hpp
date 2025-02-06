//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_eq_mult_initializer_hpp__
#define __fatrop_ip_algorithm_ip_eq_mult_initializer_hpp__
#include "fatrop/ip_algorithm/fwd.hpp"
#include "fatrop/linear_algebra/vector.hpp"

#include "fatrop/nlp/fwd.hpp"
#include <memory>

namespace fatrop
{
    /**
     * @brief Base class for equality multiplier initializers in interior point algorithms.
     */
    class IpEqMultInitializerBase
    {
    public:
        /**
         * @brief Initialize the equality multipliers.
         */
        virtual void initialize_eq_mult() = 0;

        /**
         * @brief Reset the initializer to its initial state.
         */
        virtual void reset() = 0;
    protected:
        virtual ~IpEqMultInitializerBase() = default;
    };

    /**
     * @brief Concrete implementation of equality multiplier initializer for a specific problem type.
     * 
     * @tparam ProblemType The type of optimization problem being solved.
     */
    template <typename ProblemType> class IpEqMultInitializer : public IpEqMultInitializerBase
    {
        typedef std::shared_ptr<PdSolverOrig<ProblemType>> PdSolverSp;
        typedef std::shared_ptr<IpData<ProblemType>> IpDataSp;
        typedef IpIterate<ProblemType> IpIterateType;
        typedef std::shared_ptr<IpIterateType> IpIterateSp;

    public:
        /**
         * @brief Construct a new IpEqMultInitializer object.
         * 
         * @param ipdata Shared pointer to the interior point algorithm data.
         * @param linear_solver Shared pointer to the primal-dual linear solver.
         */
        IpEqMultInitializer(const IpDataSp &ipdata, const PdSolverSp &linear_solver);

        void initialize_eq_mult() override;
        void reset() override;

    private:
        IpDataSp ipdata_;          ///< Interior point algorithm data
        PdSolverSp linear_solver_; ///< Primal-dual linear solver
        VecRealAllocated rhs_x_;   ///< Right-hand side for primal variables
        VecRealAllocated rhs_s_;   ///< Right-hand side for slack variables
        VecRealAllocated rhs_g_;   ///< Right-hand side for equality constraints
        VecRealAllocated rhs_cl_;  ///< Right-hand side for lower bound constraints
        VecRealAllocated rhs_cu_;  ///< Right-hand side for upper bound constraints
        VecRealAllocated Dx_;      ///< Step in primal variables
        VecRealAllocated Ds_;      ///< Step in slack variables
        VecRealAllocated Deq_;     ///< Step in equality multipliers
        VecRealAllocated dummy_s_; ///< Dummy slack variables
        VecRealAllocated dummy_z_; ///< Dummy bound multipliers
        Scalar lam_max_ = 1e3;     ///< Maximum allowed value for multipliers
    };
} // namespace fatrop

#endif // __fatrop_ip_algorithm_ip_eq_mult_initializer_hpp__
