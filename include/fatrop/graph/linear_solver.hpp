//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_graph_linear_solver_hpp__
#define __fatrop_graph_linear_solver_hpp__

#include "fatrop/context/context.hpp"
#include "fatrop/graph/block_sparsity.hpp"
#include "fatrop/graph/elimination_order.hpp"
#include "fatrop/graph/fwd.hpp"
#include "fatrop/graph/linear_system.hpp"
#include "fatrop/graph/symbolic_factorization.hpp"
#include "fatrop/graph/type.hpp"
#include "fatrop/linear_algebra/linear_solver.hpp"
#include "fatrop/linear_algebra/linear_solver_return_flags.hpp"
#include "fatrop/linear_algebra/matrix.hpp"
#include "fatrop/linear_algebra/vector.hpp"

#include <vector>

namespace fatrop
{
    /**
     * @brief Block-supernodal Cholesky solver for symmetric positive definite
     *        block-sparse linear systems.
     *
     * The solver carries:
     *   - a (block-)minimum-degree elimination order computed once at
     *     construction time (symbolic phase 1);
     *   - the elimination tree and the column-wise fill pattern of the
     *     Cholesky factor (symbolic phase 2);
     *   - persistent storage for the lower-triangular factor @c L, laid out
     *     as one BLASFEO @c MAT per (filled) block.
     *
     * Calling @c solve_in_place / @c solve_in_place_rhs from the
     * @ref LinearSolver base class drives a numerical refactorisation and a
     * forward/backward substitution. Iterative refinement is provided by the
     * base class. @c solve_once recomputes the numerical factorisation;
     * @c solve_rhs reuses the factor.
     */
    class BlockCholeskySolver : public LinearSolver<BlockCholeskySolver, GraphType>
    {
    public:
        explicit BlockCholeskySolver(const BlockSparsityPattern &sp);

        LinsolReturnFlag solve_once_impl(LinearSystem<GraphType> &ls, VecRealView &x);
        void solve_rhs_impl(LinearSystem<GraphType> &ls, VecRealView &x);

        /// Tolerance below which a diagonal pivot is considered degenerate
        /// (returns @c LinsolReturnFlag::INDEFINITE).
        void set_pivot_tol(Scalar value) { pivot_tol_ = value; }

        const BlockEliminationOrder &order() const { return order_; }
        const BlockSymbolicFactorization &symbolic() const { return symbolic_; }

    private:
        // Permute / unpermute a flat vector (length total_size) by mapping
        // per-block segments according to the elimination order.
        // No heap allocation.
        void permute_to_order_(const VecRealView &src, VecRealView &dst) const;
        void permute_from_order_(const VecRealView &src, VecRealView &dst) const;

        // Numerical factorisation of the matrix held by @c ls into @c L_.
        // Returns SUCCESS or INDEFINITE. No heap allocation.
        LinsolReturnFlag factorize_(const BlockPdMatrix &matrix);

        // Triangular solves on @c x (in permuted block layout) using @c L_.
        // No heap allocation.
        void forward_solve_(VecRealView &x_perm);
        void backward_solve_(VecRealView &x_perm);

        // O(1) lookup of the L-block index for (i, j) in permuted coordinates
        // (i >= j). Returns -1 if the block is structurally absent.
        Index l_lookup_(Index i, Index j) const;

        const BlockSparsityPattern &sp_;
        BlockEliminationOrder order_;
        BlockSymbolicFactorization symbolic_;

        // Block sizes / offsets in permuted order. Per-block: size matches
        // sp_.block_size(perm_[k]); per-block offset is a prefix sum over the
        // permuted sizes.
        std::vector<Index> perm_block_sizes_;
        std::vector<Index> perm_block_offsets_;

        // Cholesky factor @c L, one MatRealAllocated per filled lower-triangle
        // block (i, j) with i >= j in PERMUTED coordinates. The capacity is
        // sized exactly at construction time; no reallocation happens later.
        std::vector<MatRealAllocated> L_;
        // For each column j (in permuted coordinates): list of (row i, index
        // into L_) for the stored blocks in column j, sorted by row index.
        // l_col_entries_[j] always starts with the diagonal entry (i == j).
        std::vector<std::vector<std::pair<Index, Index>>> l_col_entries_;
        // Dense N x N lookup table: l_index_table_[i * N + j] = index in L_
        // for the permuted lower-triangle block (i, j) with i >= j, or -1.
        // Built at construction time; read-only thereafter.
        std::vector<Index> l_index_table_;

        // Workspace for solve: vector of length total_size in permuted layout.
        VecRealAllocated x_perm_;
        // Workspace for the right-hand side (in non-permuted layout).
        VecRealAllocated rhs_;

        Scalar pivot_tol_ = 1e-12;
        bool factorized_ = false;
    };
} // namespace fatrop

#endif // __fatrop_graph_linear_solver_hpp__
