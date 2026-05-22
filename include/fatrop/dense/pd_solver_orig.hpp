//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_dense_pd_solver_orig_hpp__
#define __fatrop_dense_pd_solver_orig_hpp__

#include "fatrop/dense/fwd.hpp"
#include "fatrop/dense/pd_system_orig.hpp"
#include "fatrop/dense/type.hpp"
#include "fatrop/ip_algorithm/pd_solver_orig.hpp"
#include "fatrop/linear_algebra/fwd.hpp"

#include <memory>

namespace fatrop
{
    /**
     * @brief Primal-dual solver for the original dense problem.
     *
     * Same reduction trick as the OCP variant: eliminate the bound multipliers
     * @c zl,zu and the slacks @c s, leaving an augmented KKT system that is
     * solved by @c AugSystemSolver<DenseType>.
     */
    template <>
    class PdSolverOrig<DenseType>
        : public LinearSolver<PdSolverOrig<DenseType>, PdSystemType<DenseType>>
    {
    public:
        PdSolverOrig(const ProblemInfo<DenseType> &info,
                     const std::shared_ptr<AugSystemSolver<DenseType>> &aug_system_solver);
        LinsolReturnFlag solve_once_impl(LinearSystem<PdSystemType<DenseType>> &ls,
                                         VecRealView &x);
        void reduce(LinearSystem<PdSystemType<DenseType>> &ls);
        void dereduce(LinearSystem<PdSystemType<DenseType>> &ls, VecRealView &x);
        void solve_rhs_impl(LinearSystem<PdSystemType<DenseType>> &ls, VecRealView &x);

    private:
        VecRealAllocated sigma_inverse_;
        VecRealAllocated ss_;
        VecRealAllocated g_ii_;
        VecRealAllocated D_ii_;
        VecRealAllocated gg_;
        VecRealAllocated x_aug_;
        VecRealAllocated mult_aug_;
        std::shared_ptr<AugSystemSolver<DenseType>> aug_system_solver_;
    };
} // namespace fatrop

#endif // __fatrop_dense_pd_solver_orig_hpp__
