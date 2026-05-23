//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_graph_aug_system_solver_hpp__
#define __fatrop_graph_aug_system_solver_hpp__

#include "fatrop/common/options.hpp"
#include "fatrop/context/context.hpp"
#include "fatrop/graph/block_pd_matrix.hpp"
#include "fatrop/graph/block_sparsity.hpp"
#include "fatrop/graph/fwd.hpp"
#include "fatrop/graph/linear_solver.hpp"
#include "fatrop/graph/linear_system.hpp"
#include "fatrop/graph/problem_type.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/linear_algebra/linear_solver_return_flags.hpp"
#include "fatrop/linear_algebra/vector.hpp"
#include "fatrop/nlp/aug_system_solver.hpp"

#include <memory>

namespace fatrop
{
    /**
     * @class AugSystemSolver<GraphProblem>
     * @brief Augmented-system solver for graph-structured NLPs.
     *
     * Solves
     * \f[
     * \begin{bmatrix} H + D_x & A_i^T \\ A_i & -D_s \end{bmatrix}
     * \begin{bmatrix} x \\ \lambda_i \end{bmatrix}
     * = -\begin{bmatrix} f \\ g_i \end{bmatrix}
     * \f]
     * by eliminating the inequality multipliers (block-local Schur complement)
     * to obtain
     * \f[
     *   (H + \sum_k A_{i,k}^T D_{s,k}^{-1} A_{i,k} + \mathrm{diag}(D_x)) \, x
     *   = -(f + \sum_k A_{i,k}^T D_{s,k}^{-1} g_{i,k}),
     * \f]
     * which is symmetric positive definite under the standard interior-point
     * assumptions and shares the user-specified block sparsity of @c H. The
     * resulting block-sparse SPD system is factorised by the
     * @ref BlockCholeskySolver. The inequality multipliers are recovered by
     * back-substitution.
     */
    template <> class AugSystemSolver<GraphProblem>
    {
    public:
        AugSystemSolver(const ProblemInfo<GraphProblem> &info);

        /// Solve the augmented system. The @c eq_mult slot is unused (graph
        /// problems have no equality multipliers); only the inequality
        /// multipliers are written.
        LinsolReturnFlag solve(const ProblemInfo<GraphProblem> &info,
                               Jacobian<GraphProblem> &jacobian, Hessian<GraphProblem> &hessian,
                               const VecRealView &D_x, const VecRealView &D_s,
                               const VecRealView &f, const VecRealView &g, VecRealView &x,
                               VecRealView &eq_mult);

        /// Solve the augmented system with an additional regularisation
        /// @c D_eq on the equality block. For graph problems @c D_eq has size
        /// 0; the call is provided to preserve symmetry with the other
        /// problem-type variants.
        LinsolReturnFlag solve(const ProblemInfo<GraphProblem> &info,
                               Jacobian<GraphProblem> &jacobian, Hessian<GraphProblem> &hessian,
                               const VecRealView &D_x, const VecRealView &D_eq,
                               const VecRealView &D_s, const VecRealView &f,
                               const VecRealView &g, VecRealView &x, VecRealView &eq_mult);

        /// Reuse the previous factorisation on a new right-hand side.
        LinsolReturnFlag solve_rhs(const ProblemInfo<GraphProblem> &info,
                                   const Jacobian<GraphProblem> &jacobian,
                                   const Hessian<GraphProblem> &hessian, const VecRealView &D_s,
                                   const VecRealView &f, const VecRealView &g, VecRealView &x,
                                   VecRealView &eq_mult);

        LinsolReturnFlag solve_rhs(const ProblemInfo<GraphProblem> &info,
                                   const Jacobian<GraphProblem> &jacobian,
                                   const Hessian<GraphProblem> &hessian,
                                   const VecRealView &D_eq, const VecRealView &D_s,
                                   const VecRealView &f, const VecRealView &g, VecRealView &x,
                                   VecRealView &eq_mult);

        void register_options(OptionRegistry &registry);

        void set_pivot_tol(const Scalar &value);

    private:
        const ProblemInfo<GraphProblem> &info_;

        // Assembled SPD block matrix M = H + sum_k A_k^T D_s_k^{-1} A_k + D_x.
        // Stored in the user-supplied block sparsity pattern.
        BlockPdMatrix M_;

        // Right-hand side of the reduced SPD system, length nx.
        VecRealAllocated rhs_;

        // Block-Cholesky solver for M.
        std::unique_ptr<BlockCholeskySolver> chol_;
    };
} // namespace fatrop

#endif // __fatrop_graph_aug_system_solver_hpp__
