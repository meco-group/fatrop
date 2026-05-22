//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//
#ifndef __fatrop_dense_pd_system_orig_hpp__
#define __fatrop_dense_pd_system_orig_hpp__
#include "fatrop/context/context.hpp"
#include "fatrop/dense/fwd.hpp"
#include "fatrop/dense/type.hpp"
#include "fatrop/ip_algorithm/pd_system_orig.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/linear_algebra/vector.hpp"

// Primal-Dual System (original) for dense problems.
//
// Mirrors LinearSystem<PdSystemType<OcpType>> exactly, but specialized on
// DenseType so the type system can pick up the dense-tagged Jacobian, Hessian
// and ProblemInfo without ambiguity.
//
// System:
//    [ H + D_x    0        A_e^T  A_i^T    0     0  ] [ x   ] = [ -f_x ]
//    [    0     D_x          0    -I      -I     I  ] [ s   ] = [ -f_s ]
//    [ A_e       0        -D_e     0       0     0  ] [ λ_e ] = [ -g_e ]
//    [ A_i      -I          0    -D_e      0     0  ] [ λ_i ] = [ -g_i ]
//    [   0     Zl_i         0      0     Sl_i    0  ] [ zl  ] = [ -cl  ]
//    [   0    -Zu_i         0      0      0    Su_i ] [ zu  ] = [ -cu  ]

namespace fatrop
{
    template <> class LinearSystem<PdSystemType<DenseType>>
    {
        friend class PdSolverOrig<DenseType>;

    public:
        LinearSystem(const ProblemInfo<DenseType> &info, Jacobian<DenseType> &jac,
                     Hessian<DenseType> &hess, const VecRealView &D_x, bool De_is_zero,
                     const VecRealView &D_e, const VecRealView &Sl_i, const VecRealView &Su_i,
                     const VecRealView &Zl_i, const VecRealView &Zu_i, VecRealView &rhs_f_x,
                     VecRealView &rhs_f_s, VecRealView &rhs_g, VecRealView &rhs_cl,
                     VecRealView &rhs_cu);

        Index m() const { return m_; }
        static Index m(const ProblemInfo<DenseType> &info);

        void get_rhs(VecRealView &out);
        void set_rhs(const VecRealView &in);
        void apply_on_right(const VecRealView &x, Scalar alpha, const VecRealView &y,
                            VecRealView &out);

    private:
        const ProblemInfo<DenseType> &info_;
        const Index m_;
        Jacobian<DenseType> &jac_;
        Hessian<DenseType> &hess_;
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

#endif // __fatrop_dense_pd_system_orig_hpp__
