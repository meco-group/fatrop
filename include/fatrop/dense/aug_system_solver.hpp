//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_dense_aug_system_solver_hpp__
#define __fatrop_dense_aug_system_solver_hpp__

#include "fatrop/common/options.hpp"
#include "fatrop/context/context.hpp"
#include "fatrop/dense/fwd.hpp"
#include "fatrop/linear_algebra/linear_solver_return_flags.hpp"
#include "fatrop/linear_algebra/lu_factorization.hpp"
#include "fatrop/linear_algebra/matrix.hpp"
#include "fatrop/nlp/aug_system_solver.hpp"

namespace fatrop
{
    /**
     * @class AugSystemSolver<DenseType>
     * @brief Augmented-system solver for dense NLP problems.
     *
     * Solves
     * \f[
     * \begin{bmatrix}
     *     H + D_x & A_e^T &  A_i^T \\
     *     A_e     & -D_e  &  0     \\
     *     A_i     & 0     &  -D_i
     * \end{bmatrix}
     * \begin{bmatrix} x \\ \lambda_e \\ \lambda_i \end{bmatrix}
     * = -\begin{bmatrix} f \\ g_e \\ g_i \end{bmatrix}
     * \f]
     *
     * This is the K = 1, nu = 0 special case of @c AugSystemSolver<OcpType>: a
     * single Schur-elimination of the equality constraints followed by a
     * Cholesky factorization of the reduced KKT block. The two-phase code path
     * mirrors the OCP implementation so the rest of the IP algorithm sees the
     * same interface.
     */
    template <> class AugSystemSolver<DenseType>
    {
    public:
        AugSystemSolver(const ProblemInfo<DenseType> &info);

        /** Solve the augmented system without equality-constraint regularization. */
        LinsolReturnFlag solve(const ProblemInfo<DenseType> &info, Jacobian<DenseType> &jacobian,
                               Hessian<DenseType> &hessian, const VecRealView &D_x,
                               const VecRealView &D_s, const VecRealView &f, const VecRealView &g,
                               VecRealView &x, VecRealView &eq_mult);

        /** Solve the augmented system with equality-constraint regularization @c D_eq. */
        LinsolReturnFlag solve(const ProblemInfo<DenseType> &info, Jacobian<DenseType> &jacobian,
                               Hessian<DenseType> &hessian, const VecRealView &D_x,
                               const VecRealView &D_eq, const VecRealView &D_s,
                               const VecRealView &f, const VecRealView &g, VecRealView &x,
                               VecRealView &eq_mult);

        /** Solve the system again for a new right-hand side (no regularization). */
        LinsolReturnFlag solve_rhs(const ProblemInfo<DenseType> &info,
                                   const Jacobian<DenseType> &jacobian,
                                   const Hessian<DenseType> &hessian, const VecRealView &D_s,
                                   const VecRealView &f, const VecRealView &g, VecRealView &x,
                                   VecRealView &eq_mult);

        /** Solve the system again for a new right-hand side (with @c D_eq regularization). */
        LinsolReturnFlag solve_rhs(const ProblemInfo<DenseType> &info,
                                   const Jacobian<DenseType> &jacobian,
                                   const Hessian<DenseType> &hessian, const VecRealView &D_eq,
                                   const VecRealView &D_s, const VecRealView &f,
                                   const VecRealView &g, VecRealView &x, VecRealView &eq_mult);

        void register_options(OptionRegistry &registry);

        void set_lu_fact_tol(const Scalar &value) { lu_fact_tol_ = value; }

    private:
        // Persistent buffers (allocated once, reused across solves).
        MatRealAllocated Ppt_;          // (nx + 1) x nx           — reduced Hessian + gradient block
        MatRealAllocated RSQrqt_tilde_; // (nx + 1) x nx           — augmented Hessian (incl. penalties)
        MatRealAllocated Ggt_stripe_;   // (nx + 1) x max(ng, 1)   — working copy of A_eq^T
        MatRealAllocated Ggt_tilde_;    // (nx + 1) x max(ng, 1)   — Schur factor of A_eq
        MatRealAllocated GgLt_;         // (nx + 1) x nx           — temporary for projection
        MatRealAllocated RSQrqt_hat_;   // (nx + 1) x nx           — projected reduced Hessian
        MatRealAllocated Llt_;          // (nx + 1) x nx           — Cholesky factor of reduced Hessian
        MatRealAllocated Ggt_ineq_;     // (nx + 1) x max(ng_ineq, 1) — scaled A_ineq^T (penalty buffer)

        // RHS-only ("solve_rhs") buffers.
        VecRealAllocated v_RSQrqt_tilde_;
        VecRealAllocated v_Ggt_stripe_;
        VecRealAllocated v_Ggt_tilde_;
        VecRealAllocated v_GgLt_;
        VecRealAllocated v_RSQrqt_hat_;
        VecRealAllocated v_Llt_;
        VecRealAllocated v_Ggt_ineq_;

        PermutationMatrix Pl_; // row permutation from the eq-Jacobian LU factorization
        PermutationMatrix Pr_; // column permutation

        Index rank_eq_ = 0;     // Effective rank of the equality Jacobian.
        Scalar lu_fact_tol_ = 1e-5;
    };
} // namespace fatrop

#endif // __fatrop_dense_aug_system_solver_hpp__
