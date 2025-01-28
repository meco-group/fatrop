//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_linear_algebra_linear_solver_hxx__
#define __fatrop_linear_algebra_linear_solver_hxx__

#include "fatrop/linear_algebra/linear_algebra.hpp"
#include "fatrop/linear_algebra/linear_system.hpp"
#include "linear_solver.hpp"
#include <limits> // for std::numeric_limits

namespace fatrop
{
    template <typename Derived, typename LsType>
    LinsolReturnFlag LinearSolver<Derived, LsType>::solve_once(LinearSystem<LsType> &ls,
                                                               VecRealView &x_in)
    {
        return derived().solve_once_impl(ls, x_in);
    }

    template <typename Derived, typename LsType>
    void LinearSolver<Derived, LsType>::solve_rhs(LinearSystem<LsType> &ls, VecRealView &x_in)
    {
        return derived().solve_rhs_impl(ls, x_in);
    }
    template <typename Derived, typename LsType>
    LinsolReturnFlag
    LinearSolver<Derived, LsType>::apply_iterative_refinement(LinearSystem<LsType> &ls)
    {
        Scalar res_prev = std::numeric_limits<Scalar>::max();
        ls.get_rhs(tmp_);
        const Scalar b_norm = norm_inf(tmp_);
        for (Index i = 0; i < max_it_ref; i++)
        {
            // compute residual residual = Ax + b
            ls.apply_on_right(x_, 1.0, tmp_, residual_);
            // compute residual inf norm
            const Scalar res_norm = norm_inf(residual_) / std::min(b_norm, 1.0);
            if (res_norm >= res_prev)
            {
                ls.set_rhs(x);
                return LinsolReturnFlag::ITREF_INCREASE;
            }
            else
            {
                // "accept" x_
                veccp(m_, x_, 0, x, 0);
            }
            if (res_norm < tol_ && i >= min_it_ref)
            {
                ls.set_rhs(x);
                return LinsolReturnFlag::SUCCESS;
            }
            if (i == max_it_ref - 1)
            {
                ls.set_rhs(x);
                return LinsolReturnFlag::ITREF_MAX_ITER;
            }
            ls.set_rhs(residual_);
            // solve
            solve_rhs(ls, x_);
            // update x
            axpy(m_, 1.0, x_, 0, x, 0, x_, 0);
            res_prev = res_norm;
        }
        ls.set_rhs(x);
        return LinsolReturnFlag::UNKNOWN;
    }

    template <typename Derived, typename LsType>
    LinsolReturnFlag LinearSolver<Derived, LsType>::solve_in_place(LinearSystem<LsType> &ls)
    {
        // reset x
        vecse(m_, 0.0, x_, 0);
        // solve
        LinsolReturnFlag ret = solve_once(ls, x_);
        if (ret != LinsolReturnFlag::SUCCESS)
        {
            return ret;
        }
        // reset residual
        vecse(m_, 0.0, residual_, 0);
        // "accept" x_
        veccp(m_, x_, 0, x, 0);
        if (max_it_ref == 0)
        {
            ls.set_rhs(x);
            return LinsolReturnFlag::SUCCESS;
        }
        return apply_iterative_refinement(ls);
    }

    template <typename Derived, typename LsType>
    LinsolReturnFlag LinearSolver<Derived, LsType>::solve_in_place_rhs(LinearSystem<LsType> &ls)
    {
        // reset x
        vecse(m_, 0.0, x_, 0);
        // solve
        solve_rhs(ls, x_);
        // reset residual
        vecse(m_, 0.0, residual_, 0);
        // "accept" x_
        veccp(m_, x_, 0, x, 0);
        if (max_it_ref == 0)
        {
            ls.set_rhs(x);
            return LinsolReturnFlag::SUCCESS;
        }
        return apply_iterative_refinement(ls);
    }
} // namespace fatrop

#endif // __fatrop_linear_algebra_linear_solver_hxx__
