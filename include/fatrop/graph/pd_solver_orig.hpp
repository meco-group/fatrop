//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_graph_pd_solver_orig_hpp__
#define __fatrop_graph_pd_solver_orig_hpp__

#include "fatrop/graph/fwd.hpp"
#include "fatrop/graph/pd_system_orig.hpp"
#include "fatrop/graph/problem_type.hpp"
#include "fatrop/ip_algorithm/pd_solver_orig.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/linear_algebra/linear_solver.hpp"

#include <memory>

namespace fatrop
{
    /**
     * @brief Primal-dual solver for the original graph-problem formulation.
     *
     * Same Schur-elimination of the bound multipliers and slacks used by the
     * dense and OCP variants; the resulting reduced augmented system has no
     * equality block and is solved by @ref AugSystemSolver<GraphProblem>
     * (which itself calls the block-Cholesky solver).
     */
    template <>
    class PdSolverOrig<GraphProblem>
        : public LinearSolver<PdSolverOrig<GraphProblem>, PdSystemType<GraphProblem>>
    {
    public:
        PdSolverOrig(const ProblemInfo<GraphProblem> &info,
                     const std::shared_ptr<AugSystemSolver<GraphProblem>> &aug_system_solver);

        LinsolReturnFlag solve_once_impl(LinearSystem<PdSystemType<GraphProblem>> &ls,
                                         VecRealView &x);
        void reduce(LinearSystem<PdSystemType<GraphProblem>> &ls);
        void dereduce(LinearSystem<PdSystemType<GraphProblem>> &ls, VecRealView &x);
        void solve_rhs_impl(LinearSystem<PdSystemType<GraphProblem>> &ls, VecRealView &x);

    private:
        VecRealAllocated sigma_inverse_;
        VecRealAllocated ss_;
        VecRealAllocated g_ii_;
        VecRealAllocated D_ii_;
        VecRealAllocated gg_;
        VecRealAllocated x_aug_;
        VecRealAllocated mult_aug_;
        std::shared_ptr<AugSystemSolver<GraphProblem>> aug_system_solver_;
    };
} // namespace fatrop

#endif // __fatrop_graph_pd_solver_orig_hpp__
