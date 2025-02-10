//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//
#ifndef __fatrop_ocp_pd_system_resto_hpp__
#define __fatrop_ocp_pd_system_resto_hpp__
#include "fatrop/context/context.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/linear_algebra/vector.hpp"
#include "fatrop/ip_algorithm/pd_system_resto.hpp"
#include "fatrop/ocp/fwd.hpp"
#include "fatrop/ocp/type.hpp"

namespace fatrop
{
    template <> class LinearSystem<PdSystemResto<OcpType>>
    {
        friend class PdSolverResto<OcpType>;

    public:
        // constructor
        LinearSystem(const ProblemInfo<OcpType> &info, Jacobian<OcpType> &jac,
                     Hessian<OcpType> &hess, const VecRealView &D_x, bool D_e_is_zero,
                     const VecRealView &D_e, const VecRealView &Sl_i, const VecRealView &Su_i,
                     const VecRealView &Zl_i, const VecRealView &Zu_i, VecRealView &rhs_f_x,
                     VecRealView &rhs_f_s, VecRealView &rhs_g, VecRealView &rhs_cl,
                     VecRealView &rhs_cu);

        /**
         * @brief Get the number of rows in the linear system.
         *
         * @return Index The number of rows.
         */
        Index m() const { return m_; };
        static Index m(const ProblemInfo<OcpType> &info);

        /**
         * @brief Get the right-hand side (RHS) of the linear system.
         *
         * @param[out] out VecRealView to store the RHS.
         */
        void get_rhs(VecRealView &out);

        /**
         * @brief Set the right-hand side (RHS) of the linear system.
         *
         * @param[in] in VecRealView containing the new RHS values.
         */
        void set_rhs(const VecRealView &in);

        /**
         * @brief Apply the system matrix A to a vector x on the right (i.e., compute Ax).
         *
         * @param[in] x VecRealView representing the input vector.
         * @param[out] out VecRealView to store the result of Ax.
         */
        void apply_on_right(const VecRealView &x, Scalar alpha, const VecRealView &y,
                            VecRealView &out);

    private:
        const ProblemInfo<OcpType> &info_;
        const Index m_;
        Jacobian<OcpType> &jac_;
        Hessian<OcpType> &hess_;
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

#endif //__fatrop_ocp_pd_system_resto_hpp__
