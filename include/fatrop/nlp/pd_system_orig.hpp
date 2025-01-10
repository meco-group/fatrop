//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_nlp_pd_system_orig__
#define __fatrop_nlp_pd_system_orig__

#include "fatrop/linear_algebra/linear_system.hpp"
#include "fatrop/nlp/fwd.hpp"

namespace fatrop
{
    template <typename ProblemType> class PdSystemType
    {
    };

    template <typename ProblemType> class LinearSystem<PdSystemType<ProblemType>>
    {
        /**
         * @brief Get the number of rows in the linear system.
         *
         * @return Index The number of rows.
         */
        Index m() const;

        /**
         * @brief Get the right-hand side (RHS) of the linear system.
         *
         * @param[out] out VecRealView to store the RHS.
         */
        void get_rhs(VecRealView &out);

        /**
         * @brief Set the right-hand side (RHS) of the linear system.
         *
         * @param[in] in VecRealView containing the new RHS values.
         */
        void set_rhs(const VecRealView &in);
        /**
         * @brief Apply the system matrix A to a vector x on the right (i.e., compute Ax).
         *
         * @param[in] x VecRealView representing the input vector.
         * @param[out] out VecRealView to store the result of Ax.
         */
        void apply_on_right(const VecRealView &x, Scalar alpha, const VecRealView& y, VecRealView &out);
    };
} // namespace fatrop

#endif //__fatrop_nlp_pd_system_orig__
