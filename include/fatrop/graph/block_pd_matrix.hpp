//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_graph_block_pd_matrix_hpp__
#define __fatrop_graph_block_pd_matrix_hpp__

#include "fatrop/context/context.hpp"
#include "fatrop/graph/block_sparsity.hpp"
#include "fatrop/graph/fwd.hpp"
#include "fatrop/linear_algebra/matrix.hpp"
#include "fatrop/linear_algebra/vector.hpp"

#include <vector>

namespace fatrop
{
    /**
     * @brief Block-sparse symmetric positive definite matrix.
     *
     * Owns one BLASFEO-backed dense block per stored entry of the
     * lower-triangular pattern of @ref BlockSparsityPattern. Only the lower
     * triangle (including diagonal) is materialised; the upper triangle is
     * implicitly the transpose. The user is expected to fill the diagonal
     * blocks with symmetric data (only the lower triangle of each diagonal
     * block is read by the solver).
     */
    class BlockPdMatrix
    {
    public:
        BlockPdMatrix(const BlockSparsityPattern &sp);

        const BlockSparsityPattern &sparsity() const { return sp_; }

        /// True if the lower-triangle block @c (i, j) with @c i >= j is
        /// structurally present.
        bool has_block(Index i, Index j) const;

        /// View of the lower-triangle block @c (i, j) with @c i >= j.
        /// Asserts that the block exists.
        MatRealView block(Index i, Index j);
        const MatRealView block(Index i, Index j) const;

        /// Zero all stored block data.
        void set_zero();

        /**
         * @brief Computes @c out = M * x + alpha * y, with @c M this symmetric
         *        block matrix.
         *
         * @c x, @c y and @c out are dense vectors of length @c total_size()
         * (concatenation of per-block segments using @c block_offset).
         */
        void apply_on_right(const VecRealView &x, Scalar alpha, const VecRealView &y,
                            VecRealView &out) const;

    private:
        const BlockSparsityPattern &sp_;
        // One block per lower-triangle structural entry (i, j) with i >= j.
        std::vector<MatRealAllocated> blocks_;
        // block_index_[i * N + j] = index in blocks_, or -1 if not stored.
        std::vector<Index> block_index_;

        Index lookup_(Index i, Index j) const;
    };
} // namespace fatrop

#endif // __fatrop_graph_block_pd_matrix_hpp__
