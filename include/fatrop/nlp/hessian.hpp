//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_nlp_hessian_hpp__
#define __fatrop_nlp_hessian_hpp__
#include "fatrop/nlp/fwd.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/context/context.hpp"

namespace fatrop
{
    template<typename ProblemType>
    struct Hessian
    {
        void apply_on_right(const ProblemInfo<ProblemType>& info, const VecRealView& x, Scalar alpha, const VecRealView& y, VecRealView& out) const;
        void get_rhs(const ProblemInfo<ProblemType>& info, VecRealView& out) const;
        void set_rhs(const ProblemInfo<ProblemType>& info, const VecRealView& in);
        void set_zero();
    };
} // namespace fatrop

#endif //__fatrop_nlp_hessian_hpp__
