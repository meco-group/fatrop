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
                             work, perm_block_offsets_[kp]);
    }
    // Traverse each supernode once; TRSM and GEMM operate on every RHS column.
    for (Index s = 0; s < symbolic_.num_supernodes(); ++s)
    {
        const Index n = sn_t_total_[s];
        if (n == 0) continue;
        const Index off = sn_off_S_[s];
        const auto &P = L_supernode_[s];
        batched::trsm_llnn(n, P, 0, kColPad, work, off, work, off);
        for (const ExtOp &e : sn_ext_ops_[s])
            batched::gemm_nn(e.ni, n, -1.0, P, e.panel_row, kColPad,
                               work, off, 1.0, work, e.x_off_i, work, e.x_off_i);
    }
    for (Index s = symbolic_.num_supernodes() - 1; s >= 0; --s)
    {
        const Index n = sn_t_total_[s];
        if (n == 0) continue;
        const Index off = sn_off_S_[s];
        const auto &P = L_supernode_[s];
        for (const ExtOp &e : sn_ext_ops_[s])
            batched::gemm_tn(e.ni, n, -1.0, P, e.panel_row, kColPad,
                                         work, e.x_off_i, 1.0, work, off, work, off);
        batched::trsm_lltn(n, P, 0, kColPad, work, off, work, off);
    }
    for (Index kp = 0; kp < sp_.num_blocks(); ++kp)
    {
        const Index orig = order_.perm()[kp];
        batched::gecp(sp_.block_size(orig), work, perm_block_offsets_[kp], X, sp_.block_offset(orig));
    }
}
