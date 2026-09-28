// Copyright (c) 2026 Lander Vanroye, KU Leuven
#ifndef FATROP_LOW_RANK_AUG_SYSTEM_SOLVER_HPP
#define FATROP_LOW_RANK_AUG_SYSTEM_SOLVER_HPP

#include "fatrop/nlp/aug_system_solver_base.hpp"
#include "fatrop/nlp/aug_system_solver.hpp"
#include "fatrop/nlp/hessian_with_low_rank.hpp"
#include "fatrop/linear_algebra/low_rank_aug_update.hpp"

namespace fatrop
{
    /// Woodbury decorator around a native augmented-system factorization.
    /// Both equality-regularization paths and subsequent scalar RHS solves use
    /// the same cached correction. The base solver computes all update columns
    /// in one matrix-RHS solve, including their equality-multiplier components.
    template <typename P> class LowRankAugSystemSolver : public AugSystemSolverBase<P>
    {
    public:
        LowRankAugSystemSolver(const ProblemInfo<P> &info,
                               std::shared_ptr<AugSystemSolver<P>> base, Index history_capacity = 6)
            : base_(std::move(base)), update_(info.number_of_tangent_variables,
                                           info.number_of_eq_constraints, history_capacity) {}

        LinsolReturnFlag solve(
            const ProblemInfo<P> &info, Jacobian<P> &jac, Hessian<P> &hess,
            const VecRealView &D_x, const VecRealView &D_s,
            const VecRealView &f, const VecRealView &g,
            VecRealView &x, VecRealView &lambda) override
        {
            auto ret = base_->solve(info, jac, hess, D_x, D_s, f, g, x, lambda);
            if (ret != LinsolReturnFlag::SUCCESS) return ret;
            ret = update_.factorize(model(hess), info.number_of_eq_constraints,
                [&](const MatRealView &F, const MatRealView &G, MatRealView &X, MatRealView &L) {
                    return base_->solve_rhs(info, jac, hess, D_s, F, G, X, L);
                });
            if (ret == LinsolReturnFlag::SUCCESS) update_.correct(x, lambda);
            return ret;
        }

        LinsolReturnFlag solve(
            const ProblemInfo<P> &info, Jacobian<P> &jac, Hessian<P> &hess,
            const VecRealView &D_x, const VecRealView &D_eq, const VecRealView &D_s,
            const VecRealView &f, const VecRealView &g,
            VecRealView &x, VecRealView &lambda) override
        {
            auto ret = base_->solve(info, jac, hess, D_x, D_eq, D_s, f, g, x, lambda);
            if (ret != LinsolReturnFlag::SUCCESS) return ret;
            ret = update_.factorize(model(hess), info.number_of_eq_constraints,
                [&](const MatRealView &F, const MatRealView &G, MatRealView &X, MatRealView &L) {
                    return base_->solve_rhs(info, jac, hess, D_eq, D_s, F, G, X, L);
                });
            if (ret == LinsolReturnFlag::SUCCESS) update_.correct(x, lambda);
            return ret;
        }

        LinsolReturnFlag solve_rhs(
            const ProblemInfo<P> &info, const Jacobian<P> &jac, const Hessian<P> &hess,
            const VecRealView &D_s, const VecRealView &f, const VecRealView &g,
            VecRealView &x, VecRealView &lambda) override
        {
            auto ret = base_->solve_rhs(info, jac, hess, D_s, f, g, x, lambda);
            if (ret == LinsolReturnFlag::SUCCESS) update_.correct(x, lambda);
            return ret;
        }

        LinsolReturnFlag solve_rhs(
            const ProblemInfo<P> &info, const Jacobian<P> &jac, const Hessian<P> &hess,
            const VecRealView &D_eq, const VecRealView &D_s,
            const VecRealView &f, const VecRealView &g,
            VecRealView &x, VecRealView &lambda) override
        {
            auto ret = base_->solve_rhs(info, jac, hess, D_eq, D_s, f, g, x, lambda);
            if (ret == LinsolReturnFlag::SUCCESS) update_.correct(x, lambda);
            return ret;
        }

    private:
        static std::shared_ptr<const LowRankHessian> model(const Hessian<P> &hess)
        {
            auto *extended = dynamic_cast<const HessianWithLowRank<P> *>(&hess);
            return extended ? extended->low_rank : nullptr;
        }
        std::shared_ptr<AugSystemSolver<P>> base_;
        LowRankAugUpdate update_;
    };
}
#endif
