//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_initializer_hpp__
#define __fatrop_ip_algorithm_ip_initializer_hpp__
#include "fatrop/context/context.hpp"
#include "fatrop/ip_algorithm/fwd.hpp"
#include <memory>
namespace fatrop
{
    /**
     * @brief Base class for initializers in interior point algorithms.
     */
    class IpInitializerBase
    {
    public:
        /**
         * @brief Initialize the algorithm's starting point.
         */
        virtual void initialize() = 0;

        /**
         * @brief Reset the initializer to its initial state.
         */
        virtual void reset() = 0;

    protected:
        virtual ~IpInitializerBase() = default;
    };

    /**
     * @brief Concrete implementation of initializer for a specific problem type.
     * 
     * @tparam ProblemType The type of optimization problem being solved.
     */
    template <typename ProblemType> class IpInitializer : public IpInitializerBase
    {
        typedef std::shared_ptr<IpEqMultInitializerBase> IpEqMultInitializerSp;
        typedef std::shared_ptr<IpData<ProblemType>> IpDataSp;

    public:
        /**
         * @brief Construct a new IpInitializer object.
         * 
         * @param ipdata Shared pointer to the interior point algorithm data.
         * @param eq_mult_initializer Shared pointer to the equality multiplier initializer.
         */
        IpInitializer(const IpDataSp ipdata, const IpEqMultInitializerSp &eq_mult_initializer);

        void initialize() override;
        void reset() override;

    private:
        /**
         * @brief Initialize the slack variables.
         */
        void initialize_slacks();

        IpDataSp ipdata_;                      ///< Interior point algorithm data
        IpEqMultInitializerSp eq_mult_initializer_; ///< Equality multiplier initializer
        Scalar bound_push = 1e-2;              ///< Bound push parameter (kappa_1)
        Scalar bound_frac = 1e-2;              ///< Bound fraction parameter (kappa_2)
    };

} // namespace fatrop

#endif // __fatrop_ip_algorithm_ip_initializer_hpp__
