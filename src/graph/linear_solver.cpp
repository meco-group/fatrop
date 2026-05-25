//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#include "fatrop/graph/linear_solver.hpp"

#include "fatrop/common/exception.hpp"
#include "fatrop/graph/block_pd_matrix.hpp"
#include "fatrop/linear_algebra/linear_algebra.hpp"
#include "fatrop/linear_algebra/linear_solver.hxx"

#include <algorithm>
#include <cstdio>

using namespace fatrop;

// Instantiate the CRTP base template for GraphType.
template class fatrop::LinearSolver<BlockCholeskySolver, GraphType>;

namespace
{
// BLASFEO panel size (D_PS). With supernode-major storage, the source block
// row offsets in gemm_nt are NOT always multiples of 4 (a block can start
// anywhere in the panel). Internally, BLASFEO's dispatch reads from
// `pC - (bi % PS) * PS` which goes BEFORE the destination's pA when the
// destination column offset is 0 — those reads land in unmapped memory and
// segfault. Avoidance: add PS leading slack columns to every supernode panel
// and store all real data at column offsets >= PS. Then `cj - (bi % PS) >= 1`,
// keeping every BLASFEO read inside the allocated panel.
constexpr int kColPad = 4;

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
      col_to_sn_(sp.num_blocks(), -1), col_to_off_in_sn_(sp.num_blocks(), -1),
      x_perm_(sp.total_size()), rhs_(sp.total_size())
{
    const Index N = sp.num_blocks();

    // Permuted block sizes/offsets.
    for (Index kp = 0; kp < N; ++kp)
        perm_block_sizes_[kp] = sp.block_size(order_.perm()[kp]);
    perm_block_offsets_ = make_offsets(perm_block_sizes_);

    // Per-supernode dimensions and per-column lookups.
    const Index num_sn = symbolic_.num_supernodes();
    sn_t_total_.assign(num_sn, 0);
    sn_ext_total_.assign(num_sn, 0);
    sn_ext_row_off_.assign(num_sn, {});
    for (Index s = 0; s < num_sn; ++s)
    {
        const Index s_start = symbolic_.supernode_start(s);
        const Index s_end = symbolic_.supernode_start(s + 1);
        Index col_acc = 0;
        for (Index k = s_start; k < s_end; ++k)
        {
            col_to_sn_[k] = s;
            col_to_off_in_sn_[k] = col_acc;
            col_acc += perm_block_sizes_[k];
        }
        sn_t_total_[s] = col_acc;

        const auto &ext_set = symbolic_.lower_pattern(s_end - 1);
        sn_ext_row_off_[s].resize(ext_set.size(), 0);
        Index row_acc = 0;
        for (size_t a = 0; a < ext_set.size(); ++a)
        {
            sn_ext_row_off_[s][a] = row_acc;
            row_acc += perm_block_sizes_[ext_set[a]];
        }
        sn_ext_total_[s] = row_acc;
    }

    // Allocate supernode panels. The +1 row is the aug-row scratch used by
    // the factor's forward-substitution fusion; it is read/written only
    // during factorize_with_rhs_ and ignored by the standalone solves. The
    // leading kColPad slack columns guard against BLASFEO gemm_nt OOB reads
    // for misaligned source row offsets (see kColPad comment above).
    L_supernode_.reserve(num_sn);
    for (Index s = 0; s < num_sn; ++s)
    {
        const Index rows = sn_t_total_[s] + sn_ext_total_[s] + 1;
        const Index cols = kColPad + std::max<Index>(1, sn_t_total_[s]);
        L_supernode_.emplace_back(rows, cols);
    }
}

BlockCholeskySolver::LBlockLoc
BlockCholeskySolver::l_loc_(Index i, Index j) const
{
    fatrop_dbg_assert(i >= j && "L stores lower triangle only");
    LBlockLoc loc;
    loc.sn = col_to_sn_[j];
    // col_off applies to PANEL columns, hence the kColPad shift; row offsets
    // are NOT padded (they index into rows of the panel directly).
    loc.col_off = kColPad + col_to_off_in_sn_[j];
    const Index s_end = symbolic_.supernode_start(loc.sn + 1);
    if (i < s_end)
    {
        // i is internal to the owning supernode — row offset is just the
        // column offset for column i.
        loc.row_off = col_to_off_in_sn_[i];
    }
    else
    {
        // i is external — binary-search the sorted ext_set of the owning
        // supernode. By the supernodal invariant lp(j) ⊆ ext_set(sn(j)) ∪
        // internal-cols-above-j, so the entry must be found.
        const auto &ext_set = symbolic_.lower_pattern(s_end - 1);
        const auto it = std::lower_bound(ext_set.begin(), ext_set.end(), i);
        const Index a = static_cast<Index>(it - ext_set.begin());
        fatrop_dbg_assert(a < static_cast<Index>(ext_set.size()) && ext_set[a] == i &&
                          "block (i, j) not in symbolic factor");
        loc.row_off = sn_t_total_[loc.sn] + sn_ext_row_off_[loc.sn][a];
    }
    return loc;
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

LinsolReturnFlag
BlockCholeskySolver::factorize_with_rhs_(const BlockPdMatrix &matrix,
                                         const VecRealView &rhs_orig)
{
    // No heap allocation in this function. All buffers (L_supernode_,
    // x_perm_, rhs_) are sized at construction time and reused.
    const Index N = sp_.num_blocks();
    const auto &inv = order_.inv_perm();

    // Permute -rhs into x_perm_. The forward substitution is folded into
    // each supernode's partial potrf via the aug-row trick; on return,
    // x_perm_ holds y = L^{-1}(-rhs).
    permute_to_order_(rhs_orig, x_perm_);
    vecsc(sp_.total_size(), -1.0, x_perm_, 0);

    // Zero all supernode panels.
    //   - Fill-in blocks not present in the input matrix start at zero.
    //   - Strict upper of the internal diagonal block (panel cells between
    //     stored blocks) must be zero before potrf_l_mn because the kernels
    //     may read those cells.
    //   - The aug row gets overwritten by rowin below, but zeroing it once
    //     is cheap.
    for (auto &P : L_supernode_)
        P = 0.0;

    // Load permuted lower-triangle of input matrix directly into supernode
    // panels at the slot l_loc_(ip, jp).
    for (Index j_orig = 0; j_orig < N; ++j_orig)
    {
        for (Index i_orig : sp_.column_pattern(j_orig))
        {
            const Index ip = inv[i_orig];
            const Index jp = inv[j_orig];
            const MatRealView src = matrix.block(i_orig, j_orig);
            if (ip >= jp)
            {
                const LBlockLoc loc = l_loc_(ip, jp);
                MatRealAllocated &dst = L_supernode_[loc.sn];
                const Index ni = perm_block_sizes_[ip];
                const Index nj = perm_block_sizes_[jp];
                gecp(ni, nj, src, 0, 0, dst, loc.row_off, loc.col_off);
            }
            else
            {
                // (ip < jp): stored as the LOWER triangle of the original
                // matrix but lands in the UPPER triangle of the permuted
                // matrix. Its transpose belongs in L at (jp, ip).
                const LBlockLoc loc = l_loc_(jp, ip);
                MatRealAllocated &dst = L_supernode_[loc.sn];
                const Index nrows_src = sp_.block_size(i_orig);
                const Index ncols_src = sp_.block_size(j_orig);
                getr(nrows_src, ncols_src, src, 0, 0, dst, loc.row_off, loc.col_off);
            }
        }
    }

    // Right-looking SUPERNODAL block Cholesky over the supernode panels.
    // For each supernode S in elimination order:
    //   1. rowin the rhs strip into the aug row at the bottom of the panel,
    //   2. one blasfeo_dpotrf_l_mn factors the t_total x t_total internal
    //      block AND trsm-solves the (ext_total + 1) x t_total bottom strip
    //      (external L rows + aug row computing y_S = L_internal^{-1}(-rhs_S)),
    //   3. extract y_S from the aug row into x_perm_, propagate it to
    //      external rhs strips with one gemv per ext block,
    //   4. schur-update — for each pair (i, j) in ext_set with i >= j, do
    //      a syrk (i == j) or gemm (i > j) writing DIRECTLY into the owning
    //      supernode's panel via l_loc_. No gather, no scatter — the
    //      panel is the storage.
    const Index num_sn = symbolic_.num_supernodes();
    for (Index s = 0; s < num_sn; ++s)
    {
        const Index s_start = symbolic_.supernode_start(s);
        const Index s_end = symbolic_.supernode_start(s + 1);
        const Index t_total = sn_t_total_[s];
        if (t_total == 0)
            continue;
        const Index ext_total = sn_ext_total_[s];
        const auto &ext_set = symbolic_.lower_pattern(s_end - 1);
        const Index n_ext = static_cast<Index>(ext_set.size());

        MatRealAllocated &P = L_supernode_[s];
        const Index aug_row = t_total + ext_total;
        const Index off_S = perm_block_offsets_[s_start];

        // (1) Load the rhs strip into the aug row, one column block at a time.
        //     Panel column = kColPad + col_to_off_in_sn_[k] (kColPad slack at start).
        for (Index k = s_start; k < s_end; ++k)
        {
            const Index nk = perm_block_sizes_[k];
            const Index col_off = kColPad + col_to_off_in_sn_[k];
            rowin(nk, 1.0, x_perm_, perm_block_offsets_[k], P, aug_row, col_off);
        }

        // (2) Partial Cholesky on the (t_total + ext_total + 1) x t_total panel,
        //     rooted at column kColPad.
        potrf_l_mn(t_total + ext_total + 1, t_total, P, 0, kColPad, P, 0, kColPad);
        if (!check_reg(t_total, &P.mat(), 0, kColPad, pivot_tol_))
            return LinsolReturnFlag::INDEFINITE;

        // (3a) Extract y_S from the aug row into x_perm_.
        for (Index k = s_start; k < s_end; ++k)
        {
            const Index nk = perm_block_sizes_[k];
            const Index col_off = kColPad + col_to_off_in_sn_[k];
            rowex(nk, 1.0, P, aug_row, col_off, x_perm_, perm_block_offsets_[k]);
        }
        // (3b) Propagate y_S into external rhs strips: x[i] -= L[i, S] * y_S.
        for (Index a = 0; a < n_ext; ++a)
        {
            const Index i = ext_set[a];
            const Index ni = perm_block_sizes_[i];
            const Index off_i = perm_block_offsets_[i];
            gemv_n(ni, t_total, -1.0, P, t_total + sn_ext_row_off_[s][a], kColPad,
                   x_perm_, off_S, 1.0, x_perm_, off_i, x_perm_, off_i);
        }

        // (4) Trailing schur — write directly into owner supernodes' panels.
        for (Index a = 0; a < n_ext; ++a)
        {
            const Index i = ext_set[a];
            const Index ni = perm_block_sizes_[i];
            const Index row_a = t_total + sn_ext_row_off_[s][a];

            // D[i] -= L_ext[i, :] * L_ext[i, :]^T  in i's owner panel.
            const LBlockLoc loc_ii = l_loc_(i, i);
            MatRealAllocated &P_ii = L_supernode_[loc_ii.sn];
            // NOTE: blasfeo_dsyrk_ln (square) silently zeros beta*C when
            // m < k and the in-place case is used. The _mn flavor is correct
            // for all m, n, k.
            syrk_ln_mn(ni, ni, t_total, -1.0, P, row_a, kColPad, P, row_a, kColPad, 1.0,
                       P_ii, loc_ii.row_off, loc_ii.col_off,
                       P_ii, loc_ii.row_off, loc_ii.col_off);

            for (Index b = 0; b < a; ++b)
            {
                const Index j = ext_set[b];
                const Index nj = perm_block_sizes_[j];
                const Index row_b = t_total + sn_ext_row_off_[s][b];
                // L[i, j] -= L_ext[i, :] * L_ext[j, :]^T  in j's owner panel.
                const LBlockLoc loc_ij = l_loc_(i, j);
                MatRealAllocated &P_ij = L_supernode_[loc_ij.sn];
                gemm_nt(ni, nj, t_total, -1.0, P, row_a, kColPad, P, row_b, kColPad, 1.0,
                        P_ij, loc_ij.row_off, loc_ij.col_off,
                        P_ij, loc_ij.row_off, loc_ij.col_off);
            }
        }
    }

    factorized_ = true;
    return LinsolReturnFlag::SUCCESS;
}

void BlockCholeskySolver::forward_solve_(VecRealView &x_perm)
{
    // No heap allocation. Walks supernode panels: one trsv per supernode
    // on the t_total x t_total internal block (rooted at column kColPad),
    // one gemv per external block.
    const Index num_sn = symbolic_.num_supernodes();
    for (Index s = 0; s < num_sn; ++s)
    {
        const Index t_total = sn_t_total_[s];
        if (t_total == 0)
            continue;
        const Index s_start = symbolic_.supernode_start(s);
        const Index s_end = symbolic_.supernode_start(s + 1);
        const Index off_S = perm_block_offsets_[s_start];
        MatRealAllocated &P = L_supernode_[s];
        trsv_lnn(t_total, P, 0, kColPad, x_perm, off_S, x_perm, off_S);
        const auto &ext_set = symbolic_.lower_pattern(s_end - 1);
        for (size_t a = 0; a < ext_set.size(); ++a)
        {
            const Index i = ext_set[a];
            const Index ni = perm_block_sizes_[i];
            const Index off_i = perm_block_offsets_[i];
            gemv_n(ni, t_total, -1.0, P, t_total + sn_ext_row_off_[s][a], kColPad,
                   x_perm, off_S, 1.0, x_perm, off_i, x_perm, off_i);
        }
    }
}

void BlockCholeskySolver::backward_solve_(VecRealView &x_perm)
{
    // No heap allocation. Walks supernode panels in reverse.
    const Index num_sn = symbolic_.num_supernodes();
    for (Index s = num_sn - 1; s >= 0; --s)
    {
        const Index t_total = sn_t_total_[s];
        if (t_total == 0)
            continue;
        const Index s_start = symbolic_.supernode_start(s);
        const Index s_end = symbolic_.supernode_start(s + 1);
        const Index off_S = perm_block_offsets_[s_start];
        MatRealAllocated &P = L_supernode_[s];
        const auto &ext_set = symbolic_.lower_pattern(s_end - 1);
        for (size_t a = 0; a < ext_set.size(); ++a)
        {
            const Index i = ext_set[a];
            const Index ni = perm_block_sizes_[i];
            const Index off_i = perm_block_offsets_[i];
            gemv_t(ni, t_total, -1.0, P, t_total + sn_ext_row_off_[s][a], kColPad,
                   x_perm, off_i, 1.0, x_perm, off_S, x_perm, off_S);
        }
        trsv_ltn(t_total, P, 0, kColPad, x_perm, off_S, x_perm, off_S);
    }
}

LinsolReturnFlag BlockCholeskySolver::solve_once_impl(LinearSystem<GraphType> &ls, VecRealView &x)
{
    // No heap allocation. factorize_with_rhs_ runs the partial Cholesky and
    // simultaneously forward-substitutes (-rhs) through every supernode's
    // panel, so on return x_perm_ already holds y = L^{-1}(-rhs). Only the
    // back-substitution and the de-permutation remain.
    ls.get_rhs(rhs_);
    LinsolReturnFlag ret = factorize_with_rhs_(ls.matrix(), rhs_);
    if (ret != LinsolReturnFlag::SUCCESS)
        return ret;
    backward_solve_(x_perm_);
    permute_from_order_(x_perm_, x);
    return LinsolReturnFlag::SUCCESS;
}

void BlockCholeskySolver::solve_rhs_impl(LinearSystem<GraphType> &ls, VecRealView &x)
{
    // No heap allocation.
    fatrop_dbg_assert(factorized_ && "solve_rhs called before a successful factorisation");
    ls.get_rhs(rhs_);
    permute_to_order_(rhs_, x_perm_);
    vecsc(sp_.total_size(), -1.0, x_perm_, 0);
    forward_solve_(x_perm_);
    backward_solve_(x_perm_);
    permute_from_order_(x_perm_, x);
}
