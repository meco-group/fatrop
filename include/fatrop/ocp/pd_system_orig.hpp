//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//
#ifndef __fatrop_ocp_pd_system_orig_hpp__
#define __fatrop_ocp_pd_system_orig_hpp__
#include "fatrop/context/context.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/linear_algebra/vector.hpp"
#include "fatrop/nlp/pd_system_orig.hpp"
#include "fatrop/ocp/fwd.hpp"
#include "fatrop/ocp/type.hpp"

// Primal-Dual System (PD System)

// Original System:
// The primal-dual system can be expressed as the following set of equations:
//
//    [ H + D_x    0        A_e^T  A_d^T  A_i^T    0     0  ] [ x   ] = [ -f_x ]
//    [    0     D_x          0      0    -I      -I     I  ] [ s   ] = [ -f_s ]
//    [ A_e       0        -D_e      0     0       0     0  ] [ λ_e ] = [ -g_e ]
//    [ A_d       0          0       0     0       0     0  ] [ λ_d ] = [ -g_d ]
//    [ A_i      -I          0       0     0     -D_i    0  ] [ λ_i ] = [ -g_i ]
//    [   0     Zl_i         0       0     0     Sl_i    0  ] [ zl  ] = [ -cl  ]
//    [   0    -Zu_i         0       0     0      0    Su_i ] [ zu  ] = [ -cu  ]
//    Here, sl and su are the lower and upper bound distance, i.e. sl = s - l and su = u - l.
//    When an inequality constraint is single-sided, the corresponding Zl (or Zu) and cl (or cu) are
//    zero.

//
// Variables:
// - \(x\), \(s\): primal variables
// - \(\lambda_e\), \(\lambda_d\), \(\lambda_i\): dual variables
// - \(z\): auxiliary variables
//
// Coefficients:
// - \(H\), \(D_x\): system matrices
// - \(A_e\), \(A_d\), \(A_i\): constraint matrices
// - \(D_e\), \(Z_i\), \(S_i\): additional system matrices
// - \(f_x\), \(f_s\), \(g_e\), \(g_d\), \(g_i\), \(c\): vector coefficients

namespace fatrop
{
    template <> class LinearSystem<PdSystemType<OcpType>>
    {
        friend class PdSolverOrig<OcpType>;

    public:
        // constructor
        LinearSystem(const ProblemInfo<OcpType> &info, Jacobian<OcpType> &jac,
                     Hessian<OcpType> &hess, const VecRealView &D_x, bool inertia_e,
                     const VecRealView &D_e, const VecRealView &D_i, const VecRealView &Sl_i,
                     const VecRealView &Su_i, const VecRealView &Zl_i, const VecRealView &Zu_i,
                     VecRealView &rhs_f_x, VecRealView &rhs_f_s, VecRealView &rhs_g,
                     VecRealView &rhs_cl, VecRealView &rhs_cu);

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
        bool inertia_e_;
        const VecRealView &D_e_;
        const VecRealView &D_i_;
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

#endif //__fatrop_ocp_pd_system_orig_hpp__
