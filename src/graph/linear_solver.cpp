//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#include "fatrop/graph/linear_solver.hpp"

#include "fatrop/common/exception.hpp"
#include "fatrop/graph/block_pd_matrix.hpp"
#include "fatrop/linear_algebra/linear_algebra.hpp"
#include "fatrop/linear_algebra/linear_solver.hxx"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <new>

using namespace fatrop;

// Instantiate the CRTP base template for GraphType.
template class fatrop::LinearSolver<BlockCholeskySolver, GraphType>;

namespace
{

std::vector<Index> make_offsets(const std::vector<Index> &sizes)
{
    std::vector<Index> off(sizes.size() + 1, 0);
    for (size_t k = 0; k < sizes.size(); ++k)
        off[k + 1] = off[k] + sizes[k];
    return off;
}
} // namespace

BlockCholeskySolver::BlockCholeskySolver(const BlockSparsityPattern &sp)
    : LinearSolver<BlockCholeskySolver, GraphType>(sp.total_size()), sp_(sp), order_(sp),
      symbolic_(sp, order_), perm_block_sizes_(sp.num_blocks()),
      col_to_sn_(sp.num_blocks(), -1), col_to_off_in_sn_(sp.num_blocks(), -1),
      trsv_ws_(sp.total_size()), x_perm_(sp.total_size()), rhs_(sp.total_size()), pivot_ref_(sp.total_size())
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
    // during factorize_with_rhs_ and ignored by the standalone solves.
    // The panels and the schur workspaces share one zero-initialised arena:
    // one allocation instead of one per matrix, and calloc leaves fresh pages
    // untouched until the first factorisation. Kernels may read past the last
    // row panel of a matrix, so the arena ends with slack.
    Index max_t = 1, max_ext = 1;
    for (Index s = 0; s < num_sn; ++s)
    {
        max_t = std::max(max_t, sn_t_total_[s]);
        max_ext = std::max(max_ext, sn_ext_total_[s]);
    }
    auto panel_rows = [&](Index s) { return sn_t_total_[s] + sn_ext_total_[s] + 1; };
    auto panel_cols = [&](Index s) { return std::max<Index>(1, sn_t_total_[s]); };
    constexpr size_t align = 64;
    size_t arena_size = MatRealMap::memsize(max_ext, max_t) +
                        MatRealMap::memsize(max_ext, max_ext) +
                        MatRealMap::memsize(8, std::max(max_t, max_ext)) + align;
    for (Index s = 0; s < num_sn; ++s)
        arena_size += MatRealMap::memsize(panel_rows(s), panel_cols(s));
    L_arena_.reset(static_cast<char *>(std::calloc(arena_size, 1)));
    if (!L_arena_)
        throw std::bad_alloc();
    const std::uintptr_t base = reinterpret_cast<std::uintptr_t>(L_arena_.get());
    char *mem = L_arena_.get() + ((base + align - 1) / align * align - base);
    L_supernode_.reserve(num_sn);
    for (Index s = 0; s < num_sn; ++s)
    {
        L_supernode_.emplace_back(panel_rows(s), panel_cols(s), mem);
        mem += MatRealMap::memsize(panel_rows(s), panel_cols(s));
    }
    ext_ws_ = std::make_unique<MatRealMap>(max_ext, max_t, mem);
    mem += MatRealMap::memsize(max_ext, max_t);
    schur_ws_ = std::make_unique<MatRealMap>(max_ext, max_ext, mem);

    // ----- precomputed op lists for the numeric phase -----------------------
    const auto &inv = order_.inv_perm();

    sn_aug_row_.assign(num_sn, 0);
    sn_off_S_.assign(num_sn, 0);
    sn_ext_ops_.assign(num_sn, {});
    sn_schur_.assign(num_sn, {});

    std::vector<Index> target_row(N, -1);
    Index scattered_sn = -1;
    std::vector<size_t> open, next_open;
    for (Index s = 0; s < num_sn; ++s)
    {
        const Index s_start = symbolic_.supernode_start(s);
        const Index s_end = symbolic_.supernode_start(s + 1);
        const Index t_total = sn_t_total_[s];
        const Index ext_total = sn_ext_total_[s];
        sn_aug_row_[s] = t_total + ext_total;
        sn_off_S_[s] = perm_block_offsets_[s_start];

        const auto &ext_set = symbolic_.lower_pattern(s_end - 1);
        const Index n_ext = static_cast<Index>(ext_set.size());

        // External strips (used by (3b), forward, backward), one per run of
        // consecutive blocks.
        for (Index a = 0; a < n_ext; ++a)
        {
            if (a > 0 && ext_set[a] == ext_set[a - 1] + 1)
            {
                sn_ext_ops_[s].back().ni += perm_block_sizes_[ext_set[a]];
                continue;
            }
            ExtOp e;
            e.ni = perm_block_sizes_[ext_set[a]];
            e.x_off_i = perm_block_offsets_[ext_set[a]];
            e.panel_row = t_total + sn_ext_row_off_[s][a];
            sn_ext_ops_[s].push_back(e);
        }

        // Schur trailing update: where the blocks (a, b), b <= a, of the dense
        // update go. Blocks are grouped into rectangles that land contiguously
        // in one panel, one gead each: rows a with consecutive destination
        // rows form a run, and runs of consecutive columns with the same rows
        // form a rectangle. A diagonal run may join a rectangle that starts
        // higher; the extra rows are upper triangle in both W and the target
        // and never read.
        auto &ops = sn_schur_[s];
        open.clear();
        for (Index b = 0; b < n_ext; ++b)
        {
            // Rows a >= b land in the target supernode owning column b, at
            // the panel rows scattered in target_row.
            const Index t = col_to_sn_[ext_set[b]];
            if (t != scattered_sn)
            {
                scattered_sn = t;
                const Index t_end = symbolic_.supernode_start(t + 1);
                for (Index k = symbolic_.supernode_start(t); k < t_end; ++k)
                    target_row[k] = col_to_off_in_sn_[k];
                const auto &t_ext = symbolic_.lower_pattern(t_end - 1);
                for (size_t a = 0; a < t_ext.size(); ++a)
                    target_row[t_ext[a]] = sn_t_total_[t] + sn_ext_row_off_[t][a];
            }
            next_open.clear();
            for (Index a = b; a < n_ext;)
            {
                SchurOp run;
                run.target_sn = t;
                run.dst_row = target_row[ext_set[a]];
                run.dst_col = col_to_off_in_sn_[ext_set[b]];
                run.ni = 0;
                run.nj = perm_block_sizes_[ext_set[b]];
                run.src_row = sn_ext_row_off_[s][a];
                run.src_col = sn_ext_row_off_[s][b];
                for (; a < n_ext && target_row[ext_set[a]] == run.dst_row + run.ni; ++a)
                    run.ni += perm_block_sizes_[ext_set[a]];
                const bool diagonal = run.src_row == run.src_col;
                bool merged = false;
                for (size_t r : open)
                {
                    SchurOp &op = ops[r];
                    if (op.target_sn == run.target_sn && op.src_col + op.nj == run.src_col &&
                        op.dst_col + op.nj == run.dst_col &&
                        op.src_row + op.ni == run.src_row + run.ni &&
                        op.dst_row - op.src_row == run.dst_row - run.src_row &&
                        (op.src_row == run.src_row || (diagonal && op.src_row < run.src_row)))
                    {
                        op.nj += run.nj;
                        next_open.push_back(r);
                        merged = true;
                        break;
                    }
                }
                if (!merged)
                {
                    next_open.push_back(ops.size());
                    ops.push_back(run);
                }
            }
            open.swap(next_open);
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
            op.src_id = static_cast<Index>(load_ops_.size());
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
    loc.col_off = col_to_off_in_sn_[j];
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
        // path uses precomputed @c sn_schur_ / @c load_ops_).
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
    // resolved during the symbolic phase (load_ops_, sn_ext_ops_, sn_schur_); this function is a flat traversal
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
        const MatRealView src = matrix.block_by_id(op.src_id);
        MatRealMap &dst = L_supernode_[op.dst_sn];
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
    //      external rhs strips with one gemv per run of ext blocks,
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
        MatRealMap &P = L_supernode_[s];
        const Index aug_row = sn_aug_row_[s];
        const Index off_S = sn_off_S_[s];

        // (1) Load the rhs strip into the aug row.
        rowin(t_total, 1.0, x_perm_, off_S, P, aug_row, 0);

        // (2) Partial Cholesky on the (t_total + ext_total + 1) x t_total panel.
        // In place, so save the diagonal check_reg needs.
        for (Index i = 0; i < t_total; ++i)
            pivot_ref_(i) = P(i, i);
        potrf_l_mn(t_total + ext_total + 1, t_total, P, 0, 0, P, 0, 0);
        if (!check_reg(t_total, P, 0, 0, pivot_ref_, 0, pivot_tol_))
            return LinsolReturnFlag::INDEFINITE;

        // (3a) Extract y_S from the aug row into x_perm_.
        rowex(t_total, 1.0, P, aug_row, 0, x_perm_, off_S);
        // (3b) Propagate y_S into external rhs strips: x[i] -= L[i, S] * y_S.
        const auto &exts = sn_ext_ops_[s];
        for (const ExtOp &e : exts)
            gemv_n(e.ni, t_total, -1.0, P, e.panel_row, 0,
                   x_perm_, off_S, 1.0, x_perm_, e.x_off_i, x_perm_, e.x_off_i);

        // (4) Trailing schur: W = L_ext L_ext^T in one syrk on a copy of the
        // external strip, then scattered into the owner supernodes' panels.
        if (ext_total > 0)
        {
            MatRealMap &E = *ext_ws_;
            MatRealMap &W = *schur_ws_;
            gecp(ext_total, t_total, P, t_total, 0, E, 0, 0);
            syrk_ln_mn(ext_total, ext_total, t_total, 1.0, E, 0, 0, E, 0, 0, 0.0, W, 0, 0, W, 0, 0);
            for (const SchurOp &op : sn_schur_[s])
                gead(op.ni, op.nj, -1.0, W, op.src_row, op.src_col,
                     L_supernode_[op.target_sn], op.dst_row, op.dst_col);
        }
    }

    factorized_ = true;
    return LinsolReturnFlag::SUCCESS;
}

void BlockCholeskySolver::forward_solve_(VecRealView &x_perm)
{
    // No heap allocation. Walks supernode panels: one trsv per supernode
    // on the t_total x t_total internal block,
    // one gemv per run of external blocks — all driven by the cached sn_ext_ops_.
    const Index num_sn = symbolic_.num_supernodes();
    for (Index s = 0; s < num_sn; ++s)
    {
        const Index t_total = sn_t_total_[s];
        if (t_total == 0)
            continue;
        const Index off_S = sn_off_S_[s];
        MatRealMap &P = L_supernode_[s];
        trsv_lnn(t_total, P, 0, 0, x_perm, off_S, x_perm, off_S);
        for (const ExtOp &e : sn_ext_ops_[s])
            gemv_n(e.ni, t_total, -1.0, P, e.panel_row, 0,
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
        MatRealMap &P = L_supernode_[s];
        for (const ExtOp &e : sn_ext_ops_[s])
            gemv_t(e.ni, t_total, -1.0, P, e.panel_row, 0,
                   x_perm, e.x_off_i, 1.0, x_perm, off_S, x_perm, off_S);
        trsv_ltn(t_total, P, 0, 0, x_perm, off_S, trsv_ws_);
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
