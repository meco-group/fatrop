//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_graph_pd_solver_resto_hpp__
#define __fatrop_graph_pd_solver_resto_hpp__

#include "fatrop/graph/fwd.hpp"
#include "fatrop/graph/problem_type.hpp"
#include "fatrop/ip_algorithm/pd_solver_resto.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/linear_algebra/linear_solver.hpp"
#include "fatrop/linear_algebra/vector.hpp"

#include <memory>

namespace fatrop
{
    template <>
    class PdSolverResto<GraphProblem>
        : public LinearSolver<PdSolverResto<GraphProblem>, PdSystemResto<GraphProblem>>
    {
    public:
        PdSolverResto(const ProblemInfo<GraphProblem> &info,
                      const std::shared_ptr<PdSolverOrig<GraphProblem>> &orig_solver);
        LinsolReturnFlag solve_once_impl(LinearSystem<PdSystemResto<GraphProblem>> &ls,
                                         VecRealView &x);
        void reduce(LinearSystem<PdSystemResto<GraphProblem>> &ls);
        void dereduce(LinearSystem<PdSystemResto<GraphProblem>> &ls, VecRealView &x);
        void solve_rhs_impl(LinearSystem<PdSystemResto<GraphProblem>> &ls, VecRealView &x);

        const std::shared_ptr<PdSolverOrig<GraphProblem>> orig_solver_;
        VecRealAllocated D_e_orig_;
        VecRealAllocated rhs_g_orig_;
        VecRealAllocated f_pp_;
        VecRealAllocated f_nn_;
        VecRealAllocated Xpm1_;
        VecRealAllocated Xnm1_;
        VecRealAllocated x_orig_;
    };
} // namespace fatrop

#endif // __fatrop_graph_pd_solver_resto_hpp__
