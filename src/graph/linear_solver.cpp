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

    // ----- precomputed op lists for the numeric phase -----------------------
    const auto &inv = order_.inv_perm();

    sn_aug_row_.assign(num_sn, 0);
    sn_off_S_.assign(num_sn, 0);
    sn_rhs_cols_.assign(num_sn, {});
    sn_ext_ops_.assign(num_sn, {});
    sn_syrk_.assign(num_sn, {});
    sn_gemm_.assign(num_sn, {});

    for (Index s = 0; s < num_sn; ++s)
    {
        const Index s_start = symbolic_.supernode_start(s);
        const Index s_end = symbolic_.supernode_start(s + 1);
        const Index t_total = sn_t_total_[s];
        const Index ext_total = sn_ext_total_[s];
        sn_aug_row_[s] = t_total + ext_total;
        sn_off_S_[s] = perm_block_offsets_[s_start];

        // Per-column rhs in/out for steps (1) and (3a).
        sn_rhs_cols_[s].reserve(static_cast<size_t>(s_end - s_start));
        for (Index k = s_start; k < s_end; ++k)
        {
            RhsLoadOp r;
            r.nk = perm_block_sizes_[k];
            r.x_off = perm_block_offsets_[k];
            r.panel_col = kColPad + col_to_off_in_sn_[k];
            sn_rhs_cols_[s].push_back(r);
        }

        const auto &ext_set = symbolic_.lower_pattern(s_end - 1);
        const Index n_ext = static_cast<Index>(ext_set.size());

        // External strips (used by (3b), forward, backward).
        sn_ext_ops_[s].reserve(static_cast<size_t>(n_ext));
        for (Index a = 0; a < n_ext; ++a)
        {
            ExtOp e;
            e.ni = perm_block_sizes_[ext_set[a]];
            e.x_off_i = perm_block_offsets_[ext_set[a]];
            e.panel_row = t_total + sn_ext_row_off_[s][a];
            sn_ext_ops_[s].push_back(e);
        }

        // Schur trailing update — one syrk per a, then a flat gemm list for
        // pairs (a, b) with b < a. All target supernode indices, row/col
        // offsets, and block sizes are resolved here.
        sn_syrk_[s].reserve(static_cast<size_t>(n_ext));
        sn_gemm_[s].reserve(static_cast<size_t>(n_ext) * (static_cast<size_t>(n_ext) - 1) / 2);
        for (Index a = 0; a < n_ext; ++a)
        {
            const Index i = ext_set[a];
            const Index ni = perm_block_sizes_[i];
            const Index row_a = t_total + sn_ext_row_off_[s][a];

            const LBlockLoc loc_ii = l_loc_(i, i);
            SyrkOp sk;
            sk.target_sn = loc_ii.sn;
            sk.dst_row = loc_ii.row_off;
            sk.dst_col = loc_ii.col_off;
            sk.ni = ni;
            sk.src_row = row_a;
            sn_syrk_[s].push_back(sk);

            for (Index b = 0; b < a; ++b)
            {
                const Index j = ext_set[b];
                const Index nj = perm_block_sizes_[j];
                const Index row_b = t_total + sn_ext_row_off_[s][b];
                const LBlockLoc loc_ij = l_loc_(i, j);
                GemmOp g;
                g.target_sn = loc_ij.sn;
                g.dst_row = loc_ij.row_off;
                g.dst_col = loc_ij.col_off;
                g.ni = ni;
                g.nj = nj;
                g.src_row_a = row_a;
                g.src_row_b = row_b;
                sn_gemm_[s].push_back(g);
            }
        }
    }

    // Input matrix load ops — one entry per stored block in the original
    // sparsity. Computes the destination supernode/offsets once so the
    // numeric load loop is a flat traversal.
    size_t cap = 0;
    for (Index j_orig = 0; j_orig < N; ++j_orig)
        cap += sp_.column_pattern(j_orig).size();
    load_ops_.reserve(cap);
    for (Index j_orig = 0; j_orig < N; ++j_orig)
    {
        for (Index i_orig : sp_.column_pattern(j_orig))
        {
            const Index ip = inv[i_orig];
            const Index jp = inv[j_orig];
            LoadOp op;
            op.src_i = i_orig;
            op.src_j = j_orig;
            if (ip >= jp)
            {
                const LBlockLoc loc = l_loc_(ip, jp);
                op.dst_sn = loc.sn;
                op.dst_row = loc.row_off;
                op.dst_col = loc.col_off;
                op.ni = perm_block_sizes_[ip];
                op.nj = perm_block_sizes_[jp];
                op.transpose = false;
            }
            else
            {
                const LBlockLoc loc = l_loc_(jp, ip);
                op.dst_sn = loc.sn;
                op.dst_row = loc.row_off;
                op.dst_col = loc.col_off;
                op.ni = sp_.block_size(i_orig);
                op.nj = sp_.block_size(j_orig);
                op.transpose = true;
            }
            load_ops_.push_back(op);
        }
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
        // supernode. Called only from the symbolic phase (the numeric hot
        // path uses precomputed @c sn_syrk_ / @c sn_gemm_ / @c load_ops_).
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
    // x_perm_, rhs_) are sized at construction time and reused. All
    // destination supernode indices, panel offsets, and block sizes are
    // resolved during the symbolic phase (load_ops_, sn_rhs_cols_,
    // sn_ext_ops_, sn_syrk_, sn_gemm_); this function is a flat traversal
    // of those op lists.

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
    // panels using the cached load_ops_ table.
    for (const LoadOp &op : load_ops_)
    {
        const MatRealView src = matrix.block(op.src_i, op.src_j);
        MatRealAllocated &dst = L_supernode_[op.dst_sn];
        if (!op.transpose)
            gecp(op.ni, op.nj, src, 0, 0, dst, op.dst_row, op.dst_col);
        else
            getr(op.ni, op.nj, src, 0, 0, dst, op.dst_row, op.dst_col);
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
    //      supernode's panel using the precomputed op list.
    const Index num_sn = symbolic_.num_supernodes();
    for (Index s = 0; s < num_sn; ++s)
    {
        const Index t_total = sn_t_total_[s];
        if (t_total == 0)
            continue;
        const Index ext_total = sn_ext_total_[s];
        MatRealAllocated &P = L_supernode_[s];
        const Index aug_row = sn_aug_row_[s];
        const Index off_S = sn_off_S_[s];

        // (1) Load the rhs strip into the aug row, one column block at a time.
        const auto &rhs_cols = sn_rhs_cols_[s];
        for (const RhsLoadOp &r : rhs_cols)
            rowin(r.nk, 1.0, x_perm_, r.x_off, P, aug_row, r.panel_col);

        // (2) Partial Cholesky on the (t_total + ext_total + 1) x t_total panel,
        //     rooted at column kColPad.
        potrf_l_mn(t_total + ext_total + 1, t_total, P, 0, kColPad, P, 0, kColPad);
        if (!check_reg(t_total, &P.mat(), 0, kColPad, pivot_tol_))
            return LinsolReturnFlag::INDEFINITE;

        // (3a) Extract y_S from the aug row into x_perm_.
        for (const RhsLoadOp &r : rhs_cols)
            rowex(r.nk, 1.0, P, aug_row, r.panel_col, x_perm_, r.x_off);
        // (3b) Propagate y_S into external rhs strips: x[i] -= L[i, S] * y_S.
        const auto &exts = sn_ext_ops_[s];
        for (const ExtOp &e : exts)
            gemv_n(e.ni, t_total, -1.0, P, e.panel_row, kColPad,
                   x_perm_, off_S, 1.0, x_perm_, e.x_off_i, x_perm_, e.x_off_i);

        // (4) Trailing schur — write directly into owner supernodes' panels.
        // NOTE: blasfeo_dsyrk_ln (square) silently zeros beta*C when m < k
        // and the in-place case is used. The _mn flavor is correct for all
        // m, n, k.
        for (const SyrkOp &sk : sn_syrk_[s])
        {
            MatRealAllocated &P_ii = L_supernode_[sk.target_sn];
            syrk_ln_mn(sk.ni, sk.ni, t_total, -1.0,
                       P, sk.src_row, kColPad, P, sk.src_row, kColPad, 1.0,
                       P_ii, sk.dst_row, sk.dst_col,
                       P_ii, sk.dst_row, sk.dst_col);
        }
        for (const GemmOp &g : sn_gemm_[s])
        {
            MatRealAllocated &P_ij = L_supernode_[g.target_sn];
            gemm_nt(g.ni, g.nj, t_total, -1.0,
                    P, g.src_row_a, kColPad, P, g.src_row_b, kColPad, 1.0,
                    P_ij, g.dst_row, g.dst_col,
                    P_ij, g.dst_row, g.dst_col);
        }
    }

    factorized_ = true;
    return LinsolReturnFlag::SUCCESS;
}

void BlockCholeskySolver::forward_solve_(VecRealView &x_perm)
{
    // No heap allocation. Walks supernode panels: one trsv per supernode
    // on the t_total x t_total internal block (rooted at column kColPad),
    // one gemv per external block — all driven by the cached sn_ext_ops_.
    const Index num_sn = symbolic_.num_supernodes();
    for (Index s = 0; s < num_sn; ++s)
    {
        const Index t_total = sn_t_total_[s];
        if (t_total == 0)
            continue;
        const Index off_S = sn_off_S_[s];
        MatRealAllocated &P = L_supernode_[s];
        trsv_lnn(t_total, P, 0, kColPad, x_perm, off_S, x_perm, off_S);
        for (const ExtOp &e : sn_ext_ops_[s])
            gemv_n(e.ni, t_total, -1.0, P, e.panel_row, kColPad,
                   x_perm, off_S, 1.0, x_perm, e.x_off_i, x_perm, e.x_off_i);
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
        const Index off_S = sn_off_S_[s];
        MatRealAllocated &P = L_supernode_[s];
        for (const ExtOp &e : sn_ext_ops_[s])
            gemv_t(e.ni, t_total, -1.0, P, e.panel_row, kColPad,
                   x_perm, e.x_off_i, 1.0, x_perm, off_S, x_perm, off_S);
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
