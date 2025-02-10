//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ocp_pd_solver_resto_hpp__
#define __fatrop_ocp_pd_solver_resto_hpp__
// Primal-Dual System (PD System)

#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/ip_algorithm/pd_solver_resto.hpp"
// #include "fatrop/ip_algorithm/pd_solver_orig.hpp"
#include "fatrop/ocp/fwd.hpp"
// #include "fatrop/ocp/pd_system_resto.hpp"
#include "fatrop/ocp/type.hpp"
#include <memory>

namespace fatrop
{

    template <>
    class PdSolverResto<OcpType> : public LinearSolver<PdSolverResto<OcpType>, PdSystemResto<OcpType>>
    {
    public:
        PdSolverResto(const ProblemInfo<OcpType>& info, const std::shared_ptr<PdSolverOrig<OcpType>>& aug_system_solver);
        LinsolReturnFlag solve_once_impl(LinearSystem<PdSystemResto<OcpType>> &ls, VecRealView &x);
        void reduce(LinearSystem<PdSystemResto<OcpType>> &ls);
        void dereduce(LinearSystem<PdSystemResto<OcpType>> &ls, VecRealView &x);
        void solve_rhs_impl(LinearSystem<PdSystemResto<OcpType>> &ls, VecRealView &x);

        const std::shared_ptr<PdSolverOrig<OcpType>> orig_solver_;
        VecRealAllocated D_e_orig_;
        VecRealAllocated rhs_g_orig_;
        VecRealAllocated f_pp_;
        VecRealAllocated f_nn_;
        VecRealAllocated Xpm1_;
        VecRealAllocated Xnm1_;
        VecRealAllocated x_orig_;
    private:
    // 
    };

} // namespace fatrop

#endif //__fatrop_ocp_pd_solver_resto_hpp__
