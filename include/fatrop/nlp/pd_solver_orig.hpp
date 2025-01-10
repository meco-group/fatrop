//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_nlp_pd_solver_orig__
#define __fatrop_nlp_pd_solver_orig__

#include "fatrop/linear_algebra/linear_solver.hpp"
#include "fatrop/nlp/pd_system_orig.hpp"
#include "fatrop/nlp/fwd.hpp"

namespace fatrop
{
    template<typename ProblemType>
    class PdSolverOrig: public LinearSolver<PdSolverOrig<ProblemType>, PdSystemType<ProblemType>>
    {
        LinsolReturnFlag solve_once_impl(PdSystemType<ProblemType> &ls, VecRealView &x);
        void solve_rhs_impl(PdSystemType<ProblemType> &ls, VecRealView &x);
    };
} // namespace fatrop

#endif //__fatrop_nlp_pd_solver_orig__
