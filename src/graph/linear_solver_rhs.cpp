// Copyright (c) 2026 Lander Vanroye, KU Leuven
// Matrix-RHS supernode traversal reuses the native block Cholesky factors.
#include "fatrop/graph/linear_solver.hpp"
#include "fatrop/linear_algebra/batched_rhs.hpp"

using namespace fatrop;

void BlockCholeskySolver::solve_rhs(const MatRealView &F, MatRealView &X)
{
    fatrop_dbg_assert(factorized_ && "solve_rhs called before a successful factorization");
    fatrop_dbg_assert(F.m() == sp_.total_size() && X.m() == sp_.total_size() && F.n() == X.n());
    if (F.n() == 0) return;
    if (F.n() > rhs_capacity_)
    {
        for (Index off = 0; off < F.n(); off += rhs_capacity_)
        {
            const Index width = std::min(rhs_capacity_, F.n() - off);
            auto f = F.block(F.m(), width, 0, off);
            auto x = X.block(X.m(), width, 0, off);
            solve_rhs(f, x);
        }
        return;
    }
    auto work = x_perm_many_->block(x_perm_many_->m(), F.n(), 0, 0);
    for (Index kp = 0; kp < sp_.num_blocks(); ++kp)
    {
        const Index orig = order_.perm()[kp];
        batched::gecpsc(sp_.block_size(orig), -1.0, F, sp_.block_offset(orig),
                             work, work_block_off_[kp]);
    }
    // Traverse each supernode once; TRSM and GEMM operate on every RHS column.
    // GEMM only gets row offsets that are multiples of kPanel: the supernode rows
    // of work are, the external factor rows come from L_ext_aligned_ when the
    // panel's are not, and a misaligned external block of work goes through T.
    auto T = rhs_block_ws_->block(rhs_block_ws_->m(), F.n(), 0, 0);
    for (Index s = 0; s < symbolic_.num_supernodes(); ++s)
    {
        const Index n = sn_t_total_[s];
        if (n == 0) continue;
        const Index off = sn_work_off_[s];
        const auto &P = L_supernode_[s];
        batched::trsm_llnn(n, P, 0, 0, work, off, work, off);
        for (const ExtOp &e : sn_ext_ops_[s])
        {
            const auto &A = e.ws_row >= 0 ? L_ext_aligned_[s] : P;
            const Index ai = e.ws_row >= 0 ? e.ws_row : e.panel_row;
            if (aligned_(e.work_off_i))
            {
                batched::gemm_nn(e.ni, n, -1.0, A, ai, 0,
                                 work, off, 1.0, work, e.work_off_i, work, e.work_off_i);
            }
            else
            {
                batched::gemm_nn(e.ni, n, -1.0, A, ai, 0, work, off, 0.0, T, 0, T, 0);
                gead(e.ni, F.n(), 1.0, T, 0, 0, work, e.work_off_i, 0);
            }
        }
    }
    for (Index s = symbolic_.num_supernodes() - 1; s >= 0; --s)
    {
        const Index n = sn_t_total_[s];
        if (n == 0) continue;
        const Index off = sn_work_off_[s];
        const auto &P = L_supernode_[s];
        for (const ExtOp &e : sn_ext_ops_[s])
        {
            const auto &A = e.ws_row >= 0 ? L_ext_aligned_[s] : P;
            const Index ai = e.ws_row >= 0 ? e.ws_row : e.panel_row;
            if (aligned_(e.work_off_i))
            {
                batched::gemm_tn(e.ni, n, -1.0, A, ai, 0,
                                 work, e.work_off_i, 1.0, work, off, work, off);
            }
            else
            {
                gecp(e.ni, F.n(), work, e.work_off_i, 0, T, 0, 0);
                batched::gemm_tn(e.ni, n, -1.0, A, ai, 0, T, 0, 1.0, work, off, work, off);
            }
        }
        batched::trsm_lltn(n, P, 0, 0, work, off, work, off);
    }
    for (Index kp = 0; kp < sp_.num_blocks(); ++kp)
    {
        const Index orig = order_.perm()[kp];
        batched::gecp(sp_.block_size(orig), work, work_block_off_[kp], X, sp_.block_offset(orig));
    }
}
