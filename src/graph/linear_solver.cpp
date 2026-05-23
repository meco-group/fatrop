//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#include "fatrop/graph/linear_solver.hpp"

#include "fatrop/common/exception.hpp"
#include "fatrop/graph/block_pd_matrix.hpp"
#include "fatrop/linear_algebra/linear_algebra.hpp"
#include "fatrop/linear_algebra/linear_solver.hxx"

#include <algorithm>

using namespace fatrop;

// Instantiate the CRTP base template for GraphType.
template class fatrop::LinearSolver<BlockCholeskySolver, GraphType>;

namespace
{
// Compose the prefix-sum offsets given a vector of sizes.
std::vector<Index> make_offsets(const std::vector<Index> &sizes)
{
    std::vector<Index> off(sizes.size() + 1, 0);
    for (size_t k = 0; k < sizes.size(); ++k)
        off[k + 1] = off[k] + sizes[k];
    return off;
}

bool check_reg(const Index m, MAT *sA, const Index ai, const Index aj, Scalar tol)
{
    for (Index i = 0; i < m; i++)
    {
        if (blasfeo_matel_wrap(sA, ai + i, aj + i) < tol)
            return false;
    }
    return true;
}
} // namespace

BlockCholeskySolver::BlockCholeskySolver(const BlockSparsityPattern &sp)
    : LinearSolver<BlockCholeskySolver, GraphType>(sp.total_size()), sp_(sp), order_(sp),
      symbolic_(sp, order_), perm_block_sizes_(sp.num_blocks()),
      l_col_entries_(sp.num_blocks()), l_index_table_(sp.num_blocks() * sp.num_blocks(), -1),
      x_perm_(sp.total_size()), rhs_(sp.total_size())
{
    const Index N = sp.num_blocks();
    // Permuted block sizes/offsets.
    for (Index kp = 0; kp < N; ++kp)
        perm_block_sizes_[kp] = sp.block_size(order_.perm()[kp]);
    perm_block_offsets_ = make_offsets(perm_block_sizes_);

    // Reserve L_ exactly so push_back never reallocates: one diagonal entry
    // per column plus the symbolic fill pattern of each column.
    Index total_entries = N;
    for (Index j = 0; j < N; ++j)
        total_entries += static_cast<Index>(symbolic_.lower_pattern(j).size());
    L_.reserve(total_entries);

    // Allocate L: diagonal + filled lower-triangle blocks per column.
    for (Index j = 0; j < N; ++j)
    {
        const Index nj = perm_block_sizes_[j];
        const Index col_count = 1 + static_cast<Index>(symbolic_.lower_pattern(j).size());
        l_col_entries_[j].reserve(col_count);

        // Diagonal entry first.
        const Index diag_idx = static_cast<Index>(L_.size());
        l_col_entries_[j].emplace_back(j, diag_idx);
        l_index_table_[j * N + j] = diag_idx;
        L_.emplace_back(nj, nj);
        for (Index i : symbolic_.lower_pattern(j))
        {
            const Index ni = perm_block_sizes_[i];
            const Index off_idx = static_cast<Index>(L_.size());
            l_col_entries_[j].emplace_back(i, off_idx);
            l_index_table_[i * N + j] = off_idx;
            L_.emplace_back(ni, nj);
        }
    }
}

Index BlockCholeskySolver::l_lookup_(Index i, Index j) const
{
    fatrop_dbg_assert(i >= j && "L stores lower triangle only");
    return l_index_table_[i * sp_.num_blocks() + j];
}

void BlockCholeskySolver::permute_to_order_(const VecRealView &src, VecRealView &dst) const
{
    const Index N = sp_.num_blocks();
    for (Index kp = 0; kp < N; ++kp)
    {
        const Index k_orig = order_.perm()[kp];
        const Index n = sp_.block_size(k_orig);
        const Index src_off = sp_.block_offset(k_orig);
        const Index dst_off = perm_block_offsets_[kp];
        veccp(n, src, src_off, dst, dst_off);
    }
}

void BlockCholeskySolver::permute_from_order_(const VecRealView &src, VecRealView &dst) const
{
    const Index N = sp_.num_blocks();
    for (Index kp = 0; kp < N; ++kp)
    {
        const Index k_orig = order_.perm()[kp];
        const Index n = sp_.block_size(k_orig);
        const Index src_off = perm_block_offsets_[kp];
        const Index dst_off = sp_.block_offset(k_orig);
        veccp(n, src, src_off, dst, dst_off);
    }
}

LinsolReturnFlag BlockCholeskySolver::factorize_(const BlockPdMatrix &matrix)
{
    // No heap allocation in this function. All buffers (L_, l_index_table_,
    // workspace) are sized at construction time and reused.
    const Index N = sp_.num_blocks();
    const auto &inv = order_.inv_perm();

    // Zero L and populate with permuted lower-triangle of the input matrix.
    for (auto &b : L_)
        b = 0.0;
    // For each stored block (i_orig, j_orig) with i_orig >= j_orig in the
    // original matrix, route it to the (ip, jp) slot of L (with ip = inv[i_orig],
    // jp = inv[j_orig]). If ip < jp we transpose; if ip == jp the diagonal
    // block goes to (jp, jp); we only need to copy the lower triangle (the
    // upper is ignored by potrf_l_mn).
    for (Index j_orig = 0; j_orig < N; ++j_orig)
    {
        for (Index i_orig : sp_.column_pattern(j_orig))
        {
            const Index ip = inv[i_orig];
            const Index jp = inv[j_orig];
            const MatRealView src = matrix.block(i_orig, j_orig);
            if (ip == jp)
            {
                const Index n = perm_block_sizes_[jp];
                const Index l_idx = l_lookup_(jp, jp);
                MatRealAllocated &dst = L_[l_idx];
                gecp(n, n, src, 0, 0, dst, 0, 0);
            }
            else if (ip > jp)
            {
                const Index l_idx = l_lookup_(ip, jp);
                fatrop_dbg_assert(l_idx >= 0 && "missing slot in symbolic factor");
                MatRealAllocated &dst = L_[l_idx];
                const Index ni = perm_block_sizes_[ip];
                const Index nj = perm_block_sizes_[jp];
                gecp(ni, nj, src, 0, 0, dst, 0, 0);
            }
            else
            {
                // (ip < jp): the original stored entry is the LOWER-triangle of
                // the original matrix but maps to the UPPER triangle of the
                // permuted matrix. Its transpose lives in the lower triangle of
                // the permuted matrix at (jp, ip).
                const Index l_idx = l_lookup_(jp, ip);
                fatrop_dbg_assert(l_idx >= 0 && "missing slot in symbolic factor");
                MatRealAllocated &dst = L_[l_idx];
                const Index ni_perm = perm_block_sizes_[jp];
                const Index nj_perm = perm_block_sizes_[ip];
                const Index nrows_src = sp_.block_size(i_orig);
                const Index ncols_src = sp_.block_size(j_orig);
                fatrop_dbg_assert(ni_perm == ncols_src && nj_perm == nrows_src);
                (void)ni_perm;
                (void)nj_perm;
                // dst (ni_perm x nj_perm) := transpose(src (nrows_src x ncols_src))
                getr(nrows_src, ncols_src, src, 0, 0, dst, 0, 0);
            }
        }
    }

    // Right-looking block Cholesky.
    for (Index k = 0; k < N; ++k)
    {
        const Index nk = perm_block_sizes_[k];
        const Index Dk_idx = l_col_entries_[k][0].second;
        MatRealAllocated &Dk = L_[Dk_idx];
        if (nk == 0)
            continue;

        // Factorise the diagonal: Dk = chol(Dk).
        potrf_l_mn(nk, nk, Dk, 0, 0, Dk, 0, 0);
        if (!check_reg(nk, &Dk.mat(), 0, 0, pivot_tol_))
            return LinsolReturnFlag::INDEFINITE;

        // Triangular solve for off-diagonal entries:
        //   L[i, k] = L[i, k] * Dk^{-T}    for each i in lower_pattern(k).
        const auto &col_k = l_col_entries_[k];
        for (size_t a = 1; a < col_k.size(); ++a)
        {
            const Index i = col_k[a].first;
            const Index l_idx = col_k[a].second;
            const Index ni = perm_block_sizes_[i];
            MatRealAllocated &Lik = L_[l_idx];
            trsm_rltn(ni, nk, 1.0, Dk, 0, 0, Lik, 0, 0, Lik, 0, 0);
        }

        // Rank-nk update of the trailing matrix:
        //   D[i] -= L[i, k] * L[i, k]^T                              for each i,
        //   L[i, j] -= L[i, k] * L[j, k]^T          for each pair (j, i) in col k
        //                                            with j < i (i.e., a in col_k
        //                                            corresponds to j (= col_k[a].first),
        //                                            b corresponds to i (= col_k[b].first)).
        for (size_t a = 1; a < col_k.size(); ++a)
        {
            const Index j = col_k[a].first;
            const Index lj_idx = col_k[a].second;
            MatRealAllocated &Ljk = L_[lj_idx];
            const Index nj = perm_block_sizes_[j];

            // Diagonal update on column j.
            const Index Dj_idx = l_col_entries_[j][0].second;
            MatRealAllocated &Dj = L_[Dj_idx];
            syrk_ln(nj, nk, -1.0, Ljk, 0, 0, Ljk, 0, 0, 1.0, Dj, 0, 0, Dj, 0, 0);

            for (size_t b = a + 1; b < col_k.size(); ++b)
            {
                const Index i = col_k[b].first; // i > j
                const Index li_idx = col_k[b].second;
                MatRealAllocated &Lik = L_[li_idx];
                const Index ni = perm_block_sizes_[i];

                const Index lij_idx = l_lookup_(i, j);
                fatrop_dbg_assert(lij_idx >= 0 &&
                                  "fill-in (i, j) not allocated — symbolic factor inconsistent");
                MatRealAllocated &Lij = L_[lij_idx];
                // L[i, j] -= L[i, k] * L[j, k]^T
                gemm_nt(ni, nj, nk, -1.0, Lik, 0, 0, Ljk, 0, 0, 1.0, Lij, 0, 0, Lij, 0, 0);
            }
        }
    }

    factorized_ = true;
    return LinsolReturnFlag::SUCCESS;
}

void BlockCholeskySolver::forward_solve_(VecRealView &x_perm)
{
    // No heap allocation.
    const Index N = sp_.num_blocks();
    for (Index k = 0; k < N; ++k)
    {
        const Index nk = perm_block_sizes_[k];
        const Index off_k = perm_block_offsets_[k];
        if (nk == 0)
            continue;
        const auto &col_k = l_col_entries_[k];
        // x[k] = Dk^{-1} x[k]
        MatRealAllocated &Dk = L_[col_k[0].second];
        trsv_lnn(nk, Dk, 0, 0, x_perm, off_k, x_perm, off_k);
        // For each i in col_k below diagonal: x[i] -= L[i, k] * x[k].
        for (size_t a = 1; a < col_k.size(); ++a)
        {
            const Index i = col_k[a].first;
            MatRealAllocated &Lik = L_[col_k[a].second];
            const Index ni = perm_block_sizes_[i];
            const Index off_i = perm_block_offsets_[i];
            gemv_n(ni, nk, -1.0, Lik, 0, 0, x_perm, off_k, 1.0, x_perm, off_i, x_perm, off_i);
        }
    }
}

void BlockCholeskySolver::backward_solve_(VecRealView &x_perm)
{
    // No heap allocation.
    const Index N = sp_.num_blocks();
    for (Index k = N - 1; k >= 0; --k)
    {
        const Index nk = perm_block_sizes_[k];
        const Index off_k = perm_block_offsets_[k];
        if (nk == 0)
            continue;
        const auto &col_k = l_col_entries_[k];
        // For each i > k in col_k: x[k] -= L[i, k]^T * x[i].
        for (size_t a = 1; a < col_k.size(); ++a)
        {
            const Index i = col_k[a].first;
            MatRealAllocated &Lik = L_[col_k[a].second];
            const Index ni = perm_block_sizes_[i];
            const Index off_i = perm_block_offsets_[i];
            gemv_t(ni, nk, -1.0, Lik, 0, 0, x_perm, off_i, 1.0, x_perm, off_k, x_perm, off_k);
        }
        // x[k] = Dk^{-T} x[k]
        MatRealAllocated &Dk = L_[col_k[0].second];
        trsv_ltn(nk, Dk, 0, 0, x_perm, off_k, x_perm, off_k);
    }
}

LinsolReturnFlag BlockCholeskySolver::solve_once_impl(LinearSystem<GraphType> &ls, VecRealView &x)
{
    // No heap allocation: factorize_ and solve_rhs_impl are both alloc-free,
    // and all workspaces were allocated at construction time.
    LinsolReturnFlag ret = factorize_(ls.matrix());
    if (ret != LinsolReturnFlag::SUCCESS)
        return ret;
    solve_rhs_impl(ls, x);
    return LinsolReturnFlag::SUCCESS;
}

void BlockCholeskySolver::solve_rhs_impl(LinearSystem<GraphType> &ls, VecRealView &x)
{
    // No heap allocation.
    fatrop_dbg_assert(factorized_ && "solve_rhs called before a successful factorisation");
    // Compute permuted RHS: we want to solve M x = -rhs.  Stage 1: load -rhs
    // into the permuted workspace.
    ls.get_rhs(rhs_);
    permute_to_order_(rhs_, x_perm_);
    vecsc(sp_.total_size(), -1.0, x_perm_, 0);
    // Forward + backward triangular solves.
    forward_solve_(x_perm_);
    backward_solve_(x_perm_);
    // Map back to the original block ordering.
    permute_from_order_(x_perm_, x);
}
