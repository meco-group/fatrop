//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ocp_aug_system_solver_hpp__
#define __fatrop_ocp_aug_system_solver_hpp__
#include "fatrop/context/context.hpp"
#include "fatrop/linear_algebra/lu_factorization.hpp"
#include "fatrop/linear_algebra/matrix.hpp"
#include "fatrop/ocp/fwd.hpp"
#include <vector>

namespace fatrop
{
    enum LinsolReturnFlag
    {
        SUCCESS = 0,     // success
        INDEFINITE = 1,  // reduced Hessian is indefinite, no descent direction found
        NOFULL_RANK = 2, // Jacobian is (numerically) not full row rank
        UNKNOWN = 3      // unknown flag
    };

    /**
     * @class OcpAugSystemSolver
     * @brief Solves a system of equations for optimal control problems with augmented system
     * structure.
     *
     * This solver handles systems of the following form:
     *
     * \f[
     * \begin{bmatrix}
     *     H + D_x & A_e^T &  A_i^T \\
     *     A_e     & -D_e  &  0     \\
     *     A_i     & 0     & -D_i
     * \end{bmatrix}
     * \begin{bmatrix}
     *     x \\
     *     \lambda_e \\
     *     \lambda_i
     * \end{bmatrix} =
     * -\begin{bmatrix}
     *     f \\
     *     g_e \\
     *     g_i
     * \end{bmatrix}
     * \f]
     *
     * where:
     * - \f$ H \f$ is the Lagrangian Hessian matrix.
     * - \f$ A_e \f$ is the Jacobian matrix of equality constraints.
     * - \f$ A_i \f$ is the Jacobian matrix of inequality constraints.
     * - \f$ D_e \f$ and \f$ D_i \f$ are small diagonal regularization matrices for equality and
     * inequality constraints, respectively.
     * - \f$ x \f$ represents the primal variables.
     * - \f$ \lambda_e \f$ and \f$ \lambda_i \f$ are the Lagrange multipliers for equality and
     * inequality constraints.
     * - \f$ f, g_e, g_i \f$ are the corresponding residual vectors.
     *
     * ### Notes:
     * - When \f$ D_e = 0 \f$, the normal solve method should be used.
     * - When \f$ D_e > 0 \f$, the `solve_jac_reg` method should be employed to handle
     * regularization.
     */
    class OcpAugSystemSolver
    {
    public:
        OcpAugSystemSolver(const ProblemInfo<OcpType> &info);
        LinsolReturnFlag solve(const ProblemInfo<OcpType> &info, Jacobian<OcpType> &jacobian,
                               Hessian<OcpType> &hessian, const VecRealView &D_x,
                               const VecRealView &D_s, const VecRealView &f, const VecRealView &g,
                               VecRealView &x, VecRealView &eq_mult);
        // LinsolReturnFlag solve(const ProblemInfo<OcpType> &info, Jacobian<OcpType> &jacobian,
        //                        Hessian<OcpType> &hessian, const VecRealView &D_eq,
        //                        const VecRealView &D_s, const VecRealView &f, const VecRealView
        //                        &g, VecRealView &x, VecRealView &eq_mult);
        LinsolReturnFlag solve_rhs(const ProblemInfo<OcpType> &info,
                                   const Jacobian<OcpType> &jacobian,
                                   const Hessian<OcpType> &hessian, const VecRealView &D_s,
                                   const VecRealView &f, const VecRealView &g, VecRealView &x,
                                   VecRealView &eq_mult);

    private:
        // temporaries, pre-allocated during construction to avoid allocation during
        // optimization
        std::vector<MatRealAllocated> Ppt;
        std::vector<MatRealAllocated> Hh;
        std::vector<MatRealAllocated> AL;
        std::vector<MatRealAllocated> RSQrqt_tilde;
        std::vector<MatRealAllocated> Ggt_stripe;
        std::vector<MatRealAllocated> Ggt_tilde;
        std::vector<MatRealAllocated> GgLt;
        std::vector<MatRealAllocated> RSQrqt_hat;
        std::vector<MatRealAllocated> Llt;
        std::vector<MatRealAllocated> Llt_shift;
        std::vector<MatRealAllocated> GgIt_tilde;
        std::vector<MatRealAllocated> GgLIt;
        std::vector<MatRealAllocated> HhIt;
        std::vector<MatRealAllocated> PpIt_hat;
        std::vector<MatRealAllocated> LlIt;
        std::vector<MatRealAllocated> Ggt_ineq_temp;
        std::vector<VecRealAllocated> v_Ppt;
        std::vector<VecRealAllocated> v_Hh;
        std::vector<VecRealAllocated> v_AL;
        std::vector<VecRealAllocated> v_RSQrqt_tilde;
        std::vector<VecRealAllocated> v_Ggt_stripe;
        std::vector<VecRealAllocated> v_Ggt_tilde;
        std::vector<VecRealAllocated> v_GgLt;
        std::vector<VecRealAllocated> v_RSQrqt_hat;
        std::vector<VecRealAllocated> v_Llt;
        std::vector<VecRealAllocated> v_Llt_shift;
        std::vector<VecRealAllocated> v_GgIt_tilde;
        std::vector<VecRealAllocated> v_GgLIt;
        std::vector<VecRealAllocated> v_HhIt;
        std::vector<VecRealAllocated> v_PpIt_hat;
        std::vector<VecRealAllocated> v_LlIt;
        std::vector<VecRealAllocated> v_Ggt_ineq_temp;
        std::vector<VecRealAllocated> v_tmp;
        std::vector<PermutationMatrix> Pl;
        std::vector<PermutationMatrix> Pr;
        std::vector<PermutationMatrix> PlI;
        std::vector<PermutationMatrix> PrI;
        std::vector<Index> gamma;
        std::vector<Index> rho;
        Index rankI = 0;
        bool it_ref = true;
        bool perturbed_mode = false;
        double perturbed_mode_param = 1e-6;
        Scalar it_ref_acc = 1e-8;
        Scalar lu_fact_tol = 1e-5;
        bool diagnostic = false;
        bool increased_accuracy = true;
    };

} // namespace fatrop

#endif //__fatrop_ocp_aug_system_solver_hpp__
