//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ocp_pd__solver_orig_hpp__
#define __fatrop_ocp_pd__solver_orig_hpp__
// Primal-Dual System (PD System)

#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/ip_algorithm/pd_solver_orig.hpp"
#include "fatrop/ocp/fwd.hpp"
#include "fatrop/ocp/pd_system_orig.hpp"
#include "fatrop/ocp/type.hpp"
#include <memory>

namespace fatrop
{

    template <>
    class PdSolverOrig<OcpType> : public LinearSolver<PdSolverOrig<OcpType>, PdSystemType<OcpType>>
    {
    public:
        PdSolverOrig(const ProblemInfo<OcpType>& info, const std::shared_ptr<OcpAugSystemSolver>& aug_system_solver);
        LinsolReturnFlag solve_once_impl(LinearSystem<PdSystemType<OcpType>> &ls, VecRealView &x);
        void reduce(LinearSystem<PdSystemType<OcpType>> &ls);
        void dereduce(LinearSystem<PdSystemType<OcpType>> &ls, VecRealView &x);
        void solve_rhs_impl(LinearSystem<PdSystemType<OcpType>> &ls, VecRealView &x);

    private:
        VecRealAllocated sigma_inverse_;
        VecRealAllocated ss_;
        VecRealAllocated g_ii_;
        VecRealAllocated D_ii_;
        VecRealAllocated gg_;
        VecRealAllocated x_aug_;
        VecRealAllocated mult_aug_;
        std::shared_ptr<OcpAugSystemSolver> aug_system_solver_;
    };

} // namespace fatrop

#endif //__fatrop_ocp_pd_solver_orig_hpp__
