//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_linear_algebra_linear_solver_hpp__
#define __fatrop_linear_algebra_linear_solver_hpp__

#include "fatrop/context/context.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/linear_algebra/vector.hpp"
#include "fatrop/linear_algebra/linear_solver_return_flags.hpp"
#include <limits> // for std::numeric_limits

namespace fatrop
{
    template <typename Derived, typename LsType>
    class LinearSolver
    {
    public:
        LinearSolver(const Index m) : m_(m), x(m), x_(m_), residual_(m_), tmp_(m_) {}
        Derived &derived() { return *static_cast<Derived *>(this); }

        // solve without iterative refinement
        LinsolReturnFlag solve_once(LinearSystem<LsType> &ls, VecRealView &x_in);

        void solve_rhs(LinearSystem<LsType> &ls, VecRealView &x_in);

        // solve with iterative refinement
        LinsolReturnFlag solve_in_place(LinearSystem<LsType> &ls);
        LinsolReturnFlag solve_in_place_rhs(LinearSystem<LsType> &ls);
        LinsolReturnFlag apply_iterative_refinement(LinearSystem<LsType> &ls);

    protected:
        const Index m_;
        Index min_it_ref = 0;
        Index max_it_ref = 10;
        Scalar tol_ = 1e-8;
        VecRealAllocated x;
        VecRealAllocated x_;
        VecRealAllocated residual_;
        VecRealAllocated tmp_;
    };
} // namespace fatrop

#endif // __fatrop_linear_algebra_linear_solver_hpp__
