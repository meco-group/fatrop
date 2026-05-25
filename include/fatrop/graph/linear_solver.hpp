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
     * Storage is supernode-major: one BLASFEO @c MAT per supernode, shape
     *   (t_total + ext_total + 1) x t_total.
     * The top @c t_total x t_total block is the supernode's internal
     * lower-triangular factor; the next @c ext_total rows hold the external
     * L strip (L[i, S] for i in @c ext_set(S)); the trailing +1 row is scratch
     * used by the aug-row forward-substitution fusion during factor and is
     * untouched during back-substitution. With this layout the per-supernode
     * gather and scatter that the per-block layout required disappear — the
     * supernode panel IS the storage.
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
        struct LBlockLoc
        {
            Index sn;      ///< owning supernode (index into @c L_supernode_)
            Index row_off; ///< scalar row offset of the block within the panel
            Index col_off; ///< scalar column offset of the block within the panel
        };
        // Locate the L-block at permuted block coordinates (i, j) with i >= j
        // in the supernode-major storage. O(1) for internal blocks; O(log
        // |ext_set|) binary search for external blocks (ext_set is sorted).
        LBlockLoc l_loc_(Index i, Index j) const;

        void permute_to_order_(const VecRealView &src, VecRealView &dst) const;
        void permute_from_order_(const VecRealView &src, VecRealView &dst) const;

        // Combined numerical factorisation + forward substitution. The
        // "aug row" trick (mirroring the dense aug_system_solver — see
        // src/dense/aug_system_solver.cpp) carries the rhs strip as an extra
        // row through each supernode's partial potrf_l_mn so the forward
        // solve is absorbed into the factor kernel.
        LinsolReturnFlag factorize_with_rhs_(const BlockPdMatrix &matrix,
                                             const VecRealView &rhs_orig);

        // Supernodal triangular solves on @c x (in permuted block layout).
        // Used by @c solve_rhs_impl to reuse the factor without redoing it.
        void forward_solve_(VecRealView &x_perm);
        void backward_solve_(VecRealView &x_perm);

        const BlockSparsityPattern &sp_;
        BlockEliminationOrder order_;
        BlockSymbolicFactorization symbolic_;

        // Block sizes / offsets in permuted order.
        std::vector<Index> perm_block_sizes_;
        std::vector<Index> perm_block_offsets_;

        // The Cholesky factor, supernode-major. One MAT per supernode, sized
        // (sn_t_total_[s] + sn_ext_total_[s] + 1) x sn_t_total_[s].
        std::vector<MatRealAllocated> L_supernode_;

        // Per-supernode dimensions:
        //   sn_t_total_[s]   = sum of internal block sizes for supernode s
        //   sn_ext_total_[s] = sum of external block sizes for supernode s
        std::vector<Index> sn_t_total_;
        std::vector<Index> sn_ext_total_;
        // sn_ext_row_off_[s][a] = scalar row offset of the a-th external
        // block (indexed against lower_pattern_[s_end - 1]) within the
        // bottom (ext_total) strip of supernode s's panel. The absolute row
        // index in the panel is sn_t_total_[s] + sn_ext_row_off_[s][a].
        std::vector<std::vector<Index>> sn_ext_row_off_;

        // Per-(permuted)-column lookup:
        //   col_to_sn_[j]        = supernode index owning column j
        //   col_to_off_in_sn_[j] = scalar column offset of column j within
        //                          its owning supernode's panel (also equal
        //                          to the row offset for blocks whose row
        //                          index equals j, i.e. diagonal blocks).
        std::vector<Index> col_to_sn_;
        std::vector<Index> col_to_off_in_sn_;

        // Workspace for solve.
        VecRealAllocated x_perm_;
        VecRealAllocated rhs_;

        Scalar pivot_tol_ = 1e-12;
        bool factorized_ = false;
    };
} // namespace fatrop

#endif // __fatrop_graph_linear_solver_hpp__
