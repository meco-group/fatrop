//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_graph_pd_system_orig_hpp__
#define __fatrop_graph_pd_system_orig_hpp__

#include "fatrop/context/context.hpp"
#include "fatrop/graph/fwd.hpp"
#include "fatrop/graph/problem_type.hpp"
#include "fatrop/ip_algorithm/pd_system_orig.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/linear_algebra/vector.hpp"

// Primal-Dual System (original) for graph-structured problems.
//
// Mirrors LinearSystem<PdSystemType<DenseType>> with the only structural
// difference that the equality multiplier block has size 0. The shared search
// direction code in @c IpSearchDirImpl is written against this same interface.

namespace fatrop
{
    template <> class LinearSystem<PdSystemType<GraphProblem>>
    {
        friend class PdSolverOrig<GraphProblem>;

    public:
        LinearSystem(const ProblemInfo<GraphProblem> &info, Jacobian<GraphProblem> &jac,
                     Hessian<GraphProblem> &hess, const VecRealView &D_x, bool De_is_zero,
                     const VecRealView &D_e, const VecRealView &Sl_i, const VecRealView &Su_i,
                     const VecRealView &Zl_i, const VecRealView &Zu_i, VecRealView &rhs_f_x,
                     VecRealView &rhs_f_s, VecRealView &rhs_g, VecRealView &rhs_cl,
                     VecRealView &rhs_cu);

        Index m() const { return m_; }
        static Index m(const ProblemInfo<GraphProblem> &info);

        void get_rhs(VecRealView &out);
        void set_rhs(const VecRealView &in);
        void apply_on_right(const VecRealView &x, Scalar alpha, const VecRealView &y,
                            VecRealView &out);

    private:
        const ProblemInfo<GraphProblem> &info_;
        const Index m_;
        Jacobian<GraphProblem> &jac_;
        Hessian<GraphProblem> &hess_;
        const VecRealView &D_x_;
        bool De_is_zero_;
        const VecRealView &D_e_;
        const VecRealView &Sl_i_;
        const VecRealView &Su_i_;
        const VecRealView &Zl_i_;
        const VecRealView &Zu_i_;
        VecRealView &rhs_f_x_;
        VecRealView &rhs_f_s_;
        VecRealView &rhs_g_;
        VecRealView &rhs_cl_;
        VecRealView &rhs_cu_;
    };
} // namespace fatrop

#endif // __fatrop_graph_pd_system_orig_hpp__
