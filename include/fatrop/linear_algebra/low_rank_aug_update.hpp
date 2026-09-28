// Copyright (c) 2026 Lander Vanroye, KU Leuven
#ifndef FATROP_LOW_RANK_AUG_UPDATE_HPP
#define FATROP_LOW_RANK_AUG_UPDATE_HPP

#include "fatrop/linear_algebra/low_rank_hessian.hpp"
#include "fatrop/linear_algebra/linear_solver_return_flags.hpp"
#include <cmath>
#include <vector>

namespace fatrop
{
    /** Woodbury correction of an already factored augmented KKT system.
     * K0 contains sigma I, regularization, all Jacobians and barrier terms.
     * Solve K0 [ZV,ZU] = [V,U;0,0] as ONE matrix RHS, then factor
     * P = I + V^T ZV and N = I - U^T (ZU - ZV P^-1 V^T ZU).
     * Both small matrices must be positive definite for the correct inertia.
     * The resulting correction is reused for iterative refinement and SOC.
     */
    class LowRankAugUpdate
    {
    public:
        LowRankAugUpdate(Index n, Index ng, Index capacity)
            : rhs_x_(n, 2 * capacity), rhs_g_(ng, 2 * capacity),
              inverse_x_(n, 2 * capacity), inverse_g_(ng, 2 * capacity),
              negative_x_(n, capacity), negative_g_(ng, capacity),
              positive_(capacity * capacity), negative_(capacity * capacity), work_(capacity) {}

        template <typename SolveMany>
        LinsolReturnFlag factorize(const std::shared_ptr<const LowRankHessian> &model,
                                  Index ng, SolveMany solve_many)
        {
            model_ = model;
            if (!model_ || model_->rank() == 0) return LinsolReturnFlag::SUCCESS;
            const Index n = model_->W.m(), r = model_->rank(), m = r / 2;
            if (n != rhs_x_.m() || r > rhs_x_.n() || ng != rhs_g_.m())
                return LinsolReturnFlag::UNKNOWN;
            auto rhs_x = rhs_x_.block(n, r, 0, 0);
            auto rhs_g = rhs_g_.block(ng, r, 0, 0);
            auto inverse_x = inverse_x_.block(n, r, 0, 0);
            auto inverse_g = inverse_g_.block(ng, r, 0, 0);
            rhs_g = 0.;
            for (Index j = 0; j < r; ++j)
                for (Index i = 0; i < n; ++i) rhs_x(i, j) = -model_->W(i, j);
            auto status = solve_many(rhs_x, rhs_g, inverse_x, inverse_g);
            if (status != LinsolReturnFlag::SUCCESS) return status;

            for (Index j = 0; j < m; ++j)
                for (Index i = j; i < m; ++i)
                {
                    Scalar value = i == j ? 1. : 0.;
                    for (Index k = 0; k < n; ++k)
                        value += model_->W(k, 2 * i) * inverse_x_(k, 2 * j);
                    positive_[i * m + j] = value;
                }
            if (!cholesky(positive_, m)) return LinsolReturnFlag::INDEFINITE;

            for (Index j = 0; j < m; ++j)
            {
                for (Index i = 0; i < m; ++i)
                {
                    work_[i] = 0.;
                    for (Index k = 0; k < n; ++k)
                        work_[i] += model_->W(k, 2 * i) * inverse_x_(k, 2 * j + 1);
                }
                solve_cholesky(positive_, work_, m);
                for (Index k = 0; k < n; ++k)
                {
                    Scalar value = inverse_x_(k, 2 * j + 1);
                    for (Index i = 0; i < m; ++i) value -= inverse_x_(k, 2 * i) * work_[i];
                    negative_x_(k, j) = value;
                }
                for (Index k = 0; k < ng; ++k)
                {
                    Scalar value = inverse_g_(k, 2 * j + 1);
                    for (Index i = 0; i < m; ++i) value -= inverse_g_(k, 2 * i) * work_[i];
                    negative_g_(k, j) = value;
                }
            }
            for (Index j = 0; j < m; ++j)
                for (Index i = j; i < m; ++i)
                {
                    Scalar value = i == j ? 1. : 0.;
                    for (Index k = 0; k < n; ++k)
                        value -= model_->W(k, 2 * i + 1) * negative_x_(k, j);
                    negative_[i * m + j] = value;
                }
            return cholesky(negative_, m) ? LinsolReturnFlag::SUCCESS : LinsolReturnFlag::INDEFINITE;
        }

        void correct(VecRealView &x, VecRealView &lambda)
        {
            if (!model_ || model_->rank() == 0) return;
            const Index m = model_->rank() / 2;
            for (Index j = 0; j < m; ++j)
            {
                work_[j] = 0.;
                for (Index i = 0; i < x.m(); ++i) work_[j] += model_->W(i, 2 * j) * x(i);
            }
            solve_cholesky(positive_, work_, m);
            for (Index j = 0; j < m; ++j)
            {
                for (Index i = 0; i < x.m(); ++i) x(i) -= inverse_x_(i, 2 * j) * work_[j];
                for (Index i = 0; i < lambda.m(); ++i) lambda(i) -= inverse_g_(i, 2 * j) * work_[j];
            }
            for (Index j = 0; j < m; ++j)
            {
                work_[j] = 0.;
                for (Index i = 0; i < x.m(); ++i) work_[j] += model_->W(i, 2 * j + 1) * x(i);
            }
            solve_cholesky(negative_, work_, m);
            for (Index j = 0; j < m; ++j)
            {
                for (Index i = 0; i < x.m(); ++i) x(i) += negative_x_(i, j) * work_[j];
                for (Index i = 0; i < lambda.m(); ++i) lambda(i) += negative_g_(i, j) * work_[j];
            }
        }

    private:
        std::shared_ptr<const LowRankHessian> model_;
        MatRealAllocated rhs_x_, rhs_g_, inverse_x_, inverse_g_, negative_x_, negative_g_;
        std::vector<Scalar> positive_, negative_, work_;

        static bool cholesky(std::vector<Scalar> &a, Index n)
        {
            for (Index j = 0; j < n; ++j)
            {
                Scalar diag = a[j * n + j];
                for (Index k = 0; k < j; ++k) diag -= a[j * n + k] * a[j * n + k];
                if (!(diag > 0.) || !std::isfinite(diag)) return false;
                a[j * n + j] = std::sqrt(diag);
                for (Index i = j + 1; i < n; ++i)
                {
                    Scalar value = a[i * n + j];
                    for (Index k = 0; k < j; ++k) value -= a[i * n + k] * a[j * n + k];
                    a[i * n + j] = value / a[j * n + j];
                }
            }
            return true;
        }
        static void solve_cholesky(const std::vector<Scalar> &a, std::vector<Scalar> &b, Index n)
        {
            for (Index i = 0; i < n; ++i)
            {
                for (Index j = 0; j < i; ++j) b[i] -= a[i * n + j] * b[j];
                b[i] /= a[i * n + i];
            }
            for (Index i = n - 1; i >= 0; --i)
            {
                for (Index j = i + 1; j < n; ++j) b[i] -= a[j * n + i] * b[j];
                b[i] /= a[i * n + i];
            }
        }
    };
}
#endif
