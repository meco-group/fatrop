//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_pd_system_resto__
#define __fatrop_ip_algorithm_pd_system_resto__

#include "fatrop/linear_algebra/linear_system.hpp"
#include "fatrop/nlp/fwd.hpp"

namespace fatrop
{
    template <typename ProblemType> class PdSystemResto
    {
    };

    template <typename ProblemType> class LinearSystem<PdSystemResto<ProblemType>>;
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
