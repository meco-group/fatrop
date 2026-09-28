// Copyright (c) 2026 Lander Vanroye, KU Leuven
#ifndef FATROP_AUG_SYSTEM_SOLVER_BASE_HPP
#define FATROP_AUG_SYSTEM_SOLVER_BASE_HPP

#include "fatrop/nlp/fwd.hpp"
#include "fatrop/linear_algebra/vector.hpp"
#include "fatrop/linear_algebra/linear_solver_return_flags.hpp"

namespace fatrop
{
    /// Factorization and repeated-RHS interface used by primal-dual solvers.
    /// Implementations may wrap a structured solver without changing PD reduction.
    template <typename P> class AugSystemSolverBase
    {
    public:
        virtual ~AugSystemSolverBase() = default;
        virtual LinsolReturnFlag solve(
            const ProblemInfo<P> &, Jacobian<P> &, Hessian<P> &,
            const VecRealView &D_x, const VecRealView &D_s,
            const VecRealView &f, const VecRealView &g,
            VecRealView &x, VecRealView &lambda) = 0;
        virtual LinsolReturnFlag solve(
            const ProblemInfo<P> &, Jacobian<P> &, Hessian<P> &,
            const VecRealView &D_x, const VecRealView &D_eq, const VecRealView &D_s,
            const VecRealView &f, const VecRealView &g,
            VecRealView &x, VecRealView &lambda) = 0;
        virtual LinsolReturnFlag solve_rhs(
            const ProblemInfo<P> &, const Jacobian<P> &, const Hessian<P> &,
            const VecRealView &D_s, const VecRealView &f, const VecRealView &g,
            VecRealView &x, VecRealView &lambda) = 0;
        virtual LinsolReturnFlag solve_rhs(
            const ProblemInfo<P> &, const Jacobian<P> &, const Hessian<P> &,
            const VecRealView &D_eq, const VecRealView &D_s,
            const VecRealView &f, const VecRealView &g,
            VecRealView &x, VecRealView &lambda) = 0;
    };
}
#endif
