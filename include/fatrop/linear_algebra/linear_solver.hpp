//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_linear_algebra_linear_solver_hpp__
#define __fatrop_linear_algebra_linear_solver_hpp__
#include "fatrop/context/context.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/linear_algebra/linear_algebra.hpp"
#include "fatrop/linear_algebra/linear_solver_return_flags.hpp"
#include "fatrop/linear_algebra/linear_system.hpp"
#include <limits> // for std::numeric_limits

namespace fatrop
{
    template <typename Derived, typename LsType> class LinearSolver
    {
    public:
        LinearSolver(const Index m) : m_(m), x_(m_), residual_(m_), tmp_(m_) {}
        Derived &derived() { return *static_cast<Derived *>(this); }
        // solve without iterative refinement
        LinsolReturnFlag solve_once(LinearSystem<LsType> &ls, VecRealView &x)
        {
            return derived().solve_once_impl(ls, x);
        }
        void solve_rhs(LinearSystem<LsType> &ls, VecRealView &x)
        {
            return derived().solve_rhs_impl(ls, x);
        }
        // solve with iterative refinement
        LinsolReturnFlag solve(LinearSystem<LsType> &ls, VecRealView &x)
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
                    return LinsolReturnFlag::ITREF_INCREASE;
                }
                else
                {
                    // "accept" x_
                    veccp(m_, x_, 0, x, 0);
                }
                if (res_norm < tol_ && i >= min_it_ref)
                {
                    return LinsolReturnFlag::SUCCESS;
                }
                if (i == max_it_ref - 1)
                {
                    return LinsolReturnFlag::ITREF_MAX_ITER;
                }
                ls.set_rhs(residual_);
                // solve
                solve_rhs(ls, x_);
                // update x
                axpy(m_, 1.0, x_, 0, x, 0, x_, 0);
                res_prev = res_norm;
            }
            return LinsolReturnFlag::UNKNOWN;
        }

    protected:
        const Index m_;
        Index min_it_ref = 0;
        Index max_it_ref = 10;
        Scalar tol_ = 1e-8;
        VecRealAllocated x_;
        VecRealAllocated residual_;
        VecRealAllocated tmp_;
    };
} // namespace fatrop

#endif // __fatrop_linear_algebra_linear_solver_hpp__