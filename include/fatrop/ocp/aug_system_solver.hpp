//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ocp_aug_system_solver_hpp__
#define __fatrop_ocp_aug_system_solver_hpp__
#include "fatrop/context/context.hpp"
#include "fatrop/nlp/aug_system_solver.hpp"
#include "fatrop/linear_algebra/lu_factorization.hpp"
#include "fatrop/linear_algebra/matrix.hpp"
#include "fatrop/ocp/fwd.hpp"
#include "fatrop/linear_algebra/linear_solver_return_flags.hpp"
#include <vector>

namespace fatrop
{

    /**
     * @class AugSystemSolver<OcpType>
     * @brief Solves a system of equations for optimal control problems with augmented system
     * structure.
     *
     * This solver handles systems of the following form:
     *
     * \f[
     * \begin{bmatrix}
     *     H + D_x & A_e^T &  A_d^T &  A_i^T \\
     *     A_e     & -D_e  &  0     &  0     \\
     *     A_d     & 0     &  0     &  0     \\
     *     A_i     & 0     &  0     &  -D_i
     * \end{bmatrix}
     * \begin{bmatrix}
     *     x \\
     *     \lambda_e \\
     *     \lambda_d \\
     *     \lambda_i
     * \end{bmatrix} =
     * -\begin{bmatrix}
     *     f \\
     *     g_e \\
     *     g_d \\
     *     g_i
     * \end{bmatrix}
     * \f]
     *
     * where:
     * - \f$ H \f$ is the Lagrangian Hessian matrix.
     * - \f$ A_e \f$ is the Jacobian matrix of equality constraints.
     * - \f$ A_d \f$ is the Jacobian matrix of dynamics constraints.
     * - \f$ A_i \f$ is the Jacobian matrix of inequality constraints.
     * - \f$ D_x \f$, \f$ D_e \f$, and \f$ D_i \f$ are diagonal regularization matrices for primal variables,
     *   equality constraints, and inequality constraints, respectively.
     * - \f$ x \f$ represents the primal variables.
     * - \f$ \lambda_e \f$, \f$ \lambda_d \f$, and \f$ \lambda_i \f$ are the Lagrange multipliers for equality,
     *   dynamics, and inequality constraints, respectively.
     * - \f$ f, g_e, g_d, g_i \f$ are the corresponding residual vectors.
     *
     * The solver uses various numerical techniques, including LU factorization and iterative refinement,
     * to efficiently solve this system while handling potential numerical issues.
     */
    template<>
    class AugSystemSolver<OcpType>
    {
    public:
        /**
         * @brief Constructs an AugSystemSolver<OcpType> object.
         * @param info Problem information for the optimal control problem.
         */
        AugSystemSolver<OcpType>(const ProblemInfo<OcpType> &info);

        /**
         * @brief Solves the augmented system without path equality constraint regularization.
         * @param info Problem information.
         * @param jacobian Jacobian of the constraints.
         * @param hessian Hessian of the Lagrangian.
         * @param D_x Diagonal regularization for primal variables.
         * @param D_s Diagonal regularization for slack variables.
         * @param f Gradient of the objective function.
         * @param g Constraint residuals.
         * @param x [out] Solution vector for primal variables.
         * @param eq_mult [out] Solution vector for equality constraint multipliers.
         * @return Status flag indicating the outcome of the solve operation.
         */
        LinsolReturnFlag solve(const ProblemInfo<OcpType> &info, Jacobian<OcpType> &jacobian,
                               Hessian<OcpType> &hessian, const VecRealView &D_x,
                               const VecRealView &D_s, const VecRealView &f, const VecRealView &g,
                               VecRealView &x, VecRealView &eq_mult);

        /**
         * @brief Solves the augmented system with path equality constraint regularization.
         * @param info Problem information.
         * @param jacobian Jacobian of the constraints.
         * @param hessian Hessian of the Lagrangian.
         * @param D_x Diagonal regularization for primal variables.
         * @param D_eq Diagonal regularization for equality constraints.
         * @param D_s Diagonal regularization for slack variables.
         * @param f Gradient of the objective function.
         * @param g Constraint residuals.
         * @param x [out] Solution vector for primal variables.
         * @param eq_mult [out] Solution vector for equality constraint multipliers.
         * @return Status flag indicating the outcome of the solve operation.
         */
        LinsolReturnFlag solve(const ProblemInfo<OcpType> &info, Jacobian<OcpType> &jacobian,
                               Hessian<OcpType> &hessian, const VecRealView &D_x,
                               const VecRealView &D_eq, const VecRealView &D_s,
                               const VecRealView &f, const VecRealView &g, VecRealView &x,
                               VecRealView &eq_mult);

        /**
         * @brief Solves the system for a new right-hand side without path equality constraint regularization.
         * @param info Problem information.
         * @param jacobian Jacobian of the constraints.
         * @param hessian Hessian of the Lagrangian.
         * @param D_s Diagonal regularization for slack variables.
         * @param f Gradient of the objective function.
         * @param g Constraint residuals.
         * @param x [out] Solution vector for primal variables.
         * @param eq_mult [out] Solution vector for equality constraint multipliers.
         * @return Status flag indicating the outcome of the solve operation.
         */
        LinsolReturnFlag solve_rhs(const ProblemInfo<OcpType> &info,
                                   const Jacobian<OcpType> &jacobian,
                                   const Hessian<OcpType> &hessian, const VecRealView &D_s,
                                   const VecRealView &f, const VecRealView &g, VecRealView &x,
                                   VecRealView &eq_mult);

        /**
         * @brief Solves the system for a new right-hand side with path equality constraint regularization.
         * @param info Problem information.
         * @param jacobian Jacobian of the constraints.
         * @param hessian Hessian of the Lagrangian.
         * @param D_eq Diagonal regularization for equality constraints.
         * @param D_s Diagonal regularization for slack variables.
         * @param f Gradient of the objective function.
         * @param g Constraint residuals.
         * @param x [out] Solution vector for primal variables.
         * @param eq_mult [out] Solution vector for equality constraint multipliers.
         * @return Status flag indicating the outcome of the solve operation.
         */
        LinsolReturnFlag solve_rhs(const ProblemInfo<OcpType> &info,
                                   const Jacobian<OcpType> &jacobian,
                                   const Hessian<OcpType> &hessian, const VecRealView &D_eq,
                                   const VecRealView &D_s, const VecRealView &f,
                                   const VecRealView &g, VecRealView &x, VecRealView &eq_mult);

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
