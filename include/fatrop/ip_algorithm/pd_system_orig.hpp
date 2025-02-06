//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_pd_system_orig__
#define __fatrop_ip_algorithm_pd_system_orig__

#include "fatrop/linear_algebra/linear_system.hpp"
#include "fatrop/nlp/fwd.hpp"

namespace fatrop
{
    /**
     * @brief Represents the primal-dual system type for a specific problem.
     * 
     * This class is a placeholder for the primal-dual system type, which is used
     * in the interior point method to solve the linear system at each iteration.
     * 
     * @tparam ProblemType The type of optimization problem being solved.
     */
    template <typename ProblemType> class PdSystemType
    {
    };

    /**
     * @brief Specialization of the LinearSystem class for the PdSystemType.
     * 
     * This class represents the linear system arising from the primal-dual formulation
     * of the interior point method for a specific problem type.
     * 
     * @tparam ProblemType The type of optimization problem being solved.
     */
    template <typename ProblemType> class LinearSystem<PdSystemType<ProblemType>>;
    // {
    //     /**
    //      * @brief Get the number of rows in the linear system.
    //      *
    //      * @return Index The number of rows.
    //      */
    //     Index m() const;

    //     /**
    //      * @brief Get the right-hand side (RHS) of the linear system.
    //      *
    //      * @param[out] out VecRealView to store the RHS.
    //      */
    //     void get_rhs(VecRealView &out);

    //     /**
    //      * @brief Set the right-hand side (RHS) of the linear system.
    //      *
    //      * @param[in] in VecRealView containing the new RHS values.
    //      */
    //     void set_rhs(const VecRealView &in);

    //     /**
    //      * @brief Apply the system matrix A to a vector x on the right (i.e., compute Ax).
    //      *
    //      * @param[in] x VecRealView representing the input vector.
    //      * @param alpha Scalar multiplier for the result.
    //      * @param[in] y VecRealView representing an additional vector to be added.
    //      * @param[out] out VecRealView to store the result of alpha * Ax + y.
    //      */
    //     void apply_on_right(const VecRealView &x, Scalar alpha, const VecRealView& y, VecRealView &out);
    // };
} // namespace fatrop

#endif //__fatrop_ip_algorithm_pd_system_orig__
