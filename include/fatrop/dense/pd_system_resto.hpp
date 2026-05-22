//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//
#ifndef __fatrop_dense_pd_system_resto_hpp__
#define __fatrop_dense_pd_system_resto_hpp__
#include "fatrop/context/context.hpp"
#include "fatrop/dense/fwd.hpp"
#include "fatrop/dense/type.hpp"
#include "fatrop/ip_algorithm/pd_system_resto.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/linear_algebra/vector.hpp"

namespace fatrop
{
    template <> class LinearSystem<PdSystemResto<DenseType>>
    {
        friend class PdSolverResto<DenseType>;

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

#endif // __fatrop_dense_pd_system_resto_hpp__
