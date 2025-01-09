//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_nlp_jacobian_hpp__
#define __fatrop_nlp_jacobian_hpp__
#include "fatrop/context/context.hpp"
#include "fatrop/nlp/fwd.hpp"
#include "fatrop/linear_algebra/fwd.hpp"

namespace fatrop
{
    template<typename ProblemType>
    struct Jacobian
    {
        void apply_on_right(const ProblemInfo<ProblemType> & info, const VecRealView &x, Scalar alpha, const VecRealView& y, VecRealView &out) const;
        void transpose_apply_on_right(const ProblemInfo<ProblemType> & info, const VecRealView &mult_eq, Scalar alpha, const VecRealView& y, VecRealView &out) const;
        void get_rhs(const ProblemInfo<ProblemType> & info, VecRealView &rhs) const;
        void set_rhs(const ProblemInfo<ProblemType> & info, const VecRealView &rhs);
    };
} // namespace fatrop

#endif //__fatrop_nlp_jacobian_hpp__
