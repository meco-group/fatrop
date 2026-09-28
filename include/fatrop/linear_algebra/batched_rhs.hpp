// Copyright (c) 2026 Lander Vanroye, KU Leuven
#ifndef FATROP_LINEAR_ALGEBRA_BATCHED_RHS_HPP
#define FATROP_LINEAR_ALGEBRA_BATCHED_RHS_HPP

#include "fatrop/linear_algebra/blasfeo_operations.hpp"
#include <memory>
#include <stdexcept>
#include "fatrop/linear_algebra/linear_solver_return_flags.hpp"

namespace fatrop
{
namespace batched
{
// Matrix counterparts of the scalar RHS BLASFEO operations. Each column is one RHS;
// the RHS width is inferred from the views. Offsets count rows relative to the view.
// gemm_nn and gemm_tn retain the coefficient-matrix dimensions/order of gemv_n/t.
inline void gecp(Index m, const MatRealView &x, Index xi, MatRealView y, Index yi)
{
    fatrop_dbg_assert(x.n() == y.n());
    if (&x.mat() == &y.mat() && x.ai() + xi == y.ai() + yi && x.aj() == y.aj()) return;
    fatrop::gecp(m, x.n(), x, xi, 0, y, yi, 0);
}
inline void gecpsc(Index m, Scalar alpha, const MatRealView &x, Index xi,
                        MatRealView y, Index yi)
{
    gecp(m, x, xi, y, yi);
    gesc(m, y.n(), alpha, y, yi, 0);
}
inline void gead(Index m, Scalar alpha, const MatRealView &x, Index xi,
                const MatRealView &y, Index yi, MatRealView z, Index zi)
{
    // Call sites permit z == y, but never z == x with distinct y.
    gecp(m, y, yi, z, zi);
    fatrop::gead(m, x.n(), alpha, x, xi, 0, z, zi, 0);
}
inline void gemm_nn(Index m, Index k, Scalar alpha, const MatRealView &A,
                     Index ai, Index aj, const MatRealView &x, Index xi,
                     Scalar beta, const MatRealView &y, Index yi, MatRealView z, Index zi)
{
    fatrop_dbg_assert(x.n() == y.n() && x.n() == z.n());
    fatrop_dbg_assert(ai + m <= A.m() && aj + k <= A.n());
    fatrop_dbg_assert(xi + k <= x.m() && yi + m <= y.m() && zi + m <= z.m());
    if (m == 0 || z.n() == 0) return;
    blasfeo_dgemm_nn(m, z.n(), k, alpha, const_cast<MAT *>(&A.mat()), A.ai() + ai, A.aj() + aj,
                     const_cast<MAT *>(&x.mat()), x.ai() + xi, x.aj(), beta,
                     const_cast<MAT *>(&y.mat()), y.ai() + yi, y.aj(),
                     &z.mat(), z.ai() + zi, z.aj());
}
inline void gemm_tn(Index k, Index m, Scalar alpha, const MatRealView &A,
                               Index ai, Index aj, const MatRealView &x, Index xi,
                               Scalar beta, const MatRealView &y, Index yi,
                               MatRealView z, Index zi)
{
    fatrop_dbg_assert(x.n() == y.n() && x.n() == z.n());
    fatrop_dbg_assert(ai + k <= A.m() && aj + m <= A.n());
    fatrop_dbg_assert(xi + k <= x.m() && yi + m <= y.m() && zi + m <= z.m());
    if (m == 0 || z.n() == 0) return;
    blasfeo_dgemm_tn(m, z.n(), k, alpha, const_cast<MAT *>(&A.mat()), A.ai() + ai, A.aj() + aj,
                     const_cast<MAT *>(&x.mat()), x.ai() + xi, x.aj(), beta,
                     const_cast<MAT *>(&y.mat()), y.ai() + yi, y.aj(),
                     &z.mat(), z.ai() + zi, z.aj());
}

template <typename Kernel>
inline void trsm(Index m, const MatRealView &A, Index ai, Index aj,
                              const MatRealView &x, Index xi, MatRealView z, Index zi,
                              Kernel kernel)
{
    fatrop_dbg_assert(x.n() == z.n());
    fatrop_dbg_assert(ai + m <= A.m() && aj + m <= A.n());
    fatrop_dbg_assert(xi + m <= x.m() && zi + m <= z.m());
    if (m == 0 || z.n() == 0) return;
    kernel(m, z.n(), 1.0, const_cast<MAT *>(&A.mat()), A.ai() + ai, A.aj() + aj,
           const_cast<MAT *>(&x.mat()), x.ai() + xi, x.aj(), &z.mat(), z.ai() + zi, z.aj());
}
inline void trsm_llnn(Index m, const MatRealView &A, Index ai, Index aj,
                        const MatRealView &x, Index xi, MatRealView z, Index zi)
{
    trsm(m, A, ai, aj, x, xi, z, zi, blasfeo_dtrsm_llnn);
}
inline void trsm_lltn(Index m, const MatRealView &A, Index ai, Index aj,
                                  const MatRealView &x, Index xi, MatRealView z, Index zi)
{
    trsm(m, A, ai, aj, x, xi, z, zi, blasfeo_dtrsm_lltn);
}
inline void trsm_lutu(Index m, const MatRealView &A, Index ai, Index aj,
                                       const MatRealView &x, Index xi, MatRealView z, Index zi)
{
    trsm(m, A, ai, aj, x, xi, z, zi, blasfeo_dtrsm_lutu);
}
inline void trsm_lunu(Index m, Index n, const MatRealView &A, Index ai, Index aj,
                             const MatRealView &x, Index xi, MatRealView z, Index zi)
{
    // The rectangular upper system includes already known trailing solution rows.
    gecp(n, x, xi, z, zi);
    if (n > m)
        gemm_nn(m, n - m, -1.0, A, ai, aj + m, x, xi + m, 1.0, z, zi, z, zi);
    trsm(m, A, ai, aj, z, zi, z, zi, blasfeo_dtrsm_lunu);
}
inline void rowpe(PermutationMatrix &p, Index rank, MatRealView x, Index off)
{
    for (Index i = 0; i < rank; ++i)
        if (p[i] != i) rowsw(x.n(), x, off + i, 0, x, off + p[i], 0);
}
inline void rowpei(PermutationMatrix &p, Index rank, MatRealView x, Index off)
{
    for (Index i = rank - 1; i >= 0; --i)
        if (p[i] != i) rowsw(x.n(), x, off + i, 0, x, off + p[i], 0);
}
inline void divide_rows(Index m, const MatRealView &x, Index xi, const VecRealView &d,
                        Index di, MatRealView z, Index zi)
{
    fatrop_dbg_assert(x.n() == z.n());
    for (Index j = 0; j < z.n(); ++j)
        for (Index i = 0; i < m; ++i)
            z(zi + i, j) = x(xi + i, j) / d(di + i);
}
// Allocated exactly once by the solver constructor.
inline Index checked_capacity(Index capacity)
{
    if (capacity <= 0) throw std::invalid_argument("RHS capacity must be positive");
    return capacity;
}
inline void allocate_like(std::vector<MatRealAllocated> &work,
                          const std::vector<VecRealAllocated> &shape, Index capacity)
{
    work.reserve(shape.size());
    for (const auto &v : shape) work.emplace_back(v.m(), capacity);
}
class WorkspaceView
{
public:
    WorkspaceView(std::vector<MatRealAllocated> &storage, Index width)
        : storage_(storage), width_(width) {}
    MatRealView operator[](Index i) const
    {
        return storage_[i].block(storage_[i].m(), width_, 0, 0);
    }
private:
    std::vector<MatRealAllocated> &storage_;
    Index width_;
};
template <typename SolvePanel>
LinsolReturnFlag solve_panels(Index capacity, const MatRealView &F, const MatRealView &G,
                             MatRealView X, MatRealView Lambda, SolvePanel solve_panel)
{
    for (Index off = 0; off < F.n(); off += capacity)
    {
        const Index width = std::min(capacity, F.n() - off);
        auto f = F.block(F.m(), width, 0, off);
        auto g = G.block(G.m(), width, 0, off);
        auto x = X.block(X.m(), width, 0, off);
        auto lambda = Lambda.block(Lambda.m(), width, 0, off);
        const auto ret = solve_panel(f, g, x, lambda);
        if (ret != LinsolReturnFlag::SUCCESS) return ret;
    }
    return LinsolReturnFlag::SUCCESS;
}
inline void check_dimensions(Index nx, Index ng, const MatRealView &f, const MatRealView &g,
                              const MatRealView &x, const MatRealView &mult)
{
    fatrop_dbg_assert(f.m() == nx && x.m() == nx && g.m() == ng && mult.m() == ng);
    fatrop_dbg_assert(f.n() == g.n() && f.n() == x.n() && f.n() == mult.n());
}
} // namespace batched
} // namespace fatrop
#endif
