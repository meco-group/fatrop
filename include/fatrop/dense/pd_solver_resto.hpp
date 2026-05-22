//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_dense_pd_solver_resto_hpp__
#define __fatrop_dense_pd_solver_resto_hpp__

#include "fatrop/dense/fwd.hpp"
#include "fatrop/dense/type.hpp"
#include "fatrop/ip_algorithm/pd_solver_resto.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/linear_algebra/vector.hpp"

#include <memory>

namespace fatrop
{
    template <>
    class PdSolverResto<DenseType>
        : public LinearSolver<PdSolverResto<DenseType>, PdSystemResto<DenseType>>
    {
    public:
        PdSolverResto(const ProblemInfo<DenseType> &info,
                      const std::shared_ptr<PdSolverOrig<DenseType>> &orig_solver);
        LinsolReturnFlag solve_once_impl(LinearSystem<PdSystemResto<DenseType>> &ls,
                                         VecRealView &x);
        void reduce(LinearSystem<PdSystemResto<DenseType>> &ls);
        void dereduce(LinearSystem<PdSystemResto<DenseType>> &ls, VecRealView &x);
        void solve_rhs_impl(LinearSystem<PdSystemResto<DenseType>> &ls, VecRealView &x);

        const std::shared_ptr<PdSolverOrig<DenseType>> orig_solver_;
        VecRealAllocated D_e_orig_;
        VecRealAllocated rhs_g_orig_;
        VecRealAllocated f_pp_;
        VecRealAllocated f_nn_;
        VecRealAllocated Xpm1_;
        VecRealAllocated Xnm1_;
        VecRealAllocated x_orig_;
    };
} // namespace fatrop

#endif // __fatrop_dense_pd_solver_resto_hpp__
