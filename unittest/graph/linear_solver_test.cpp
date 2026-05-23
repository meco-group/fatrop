//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//
#include "../random_matrix.hpp"
#include "fatrop/graph/block_pd_matrix.hpp"
#include "fatrop/graph/block_sparsity.hpp"
#include "fatrop/graph/elimination_order.hpp"
#include "fatrop/graph/linear_solver.hpp"
#include "fatrop/graph/linear_system.hpp"
#include "fatrop/graph/symbolic_factorization.hpp"
#include "fatrop/linear_algebra/linear_algebra.hpp"
#include "fatrop/linear_algebra/linear_solver.hxx"

#include <algorithm>
#include <gtest/gtest.h>
#include <vector>

using namespace fatrop;

namespace
{
// Fill the lower triangle of a diagonal block with a symmetric random matrix
// shifted to be safely positive definite.
void fill_diagonal_block(MatRealView Mkk, Index n, Scalar shift)
{
    MatRealAllocated tmp = ::test::random_spd_matrix(n);
    gecp(n, n, tmp, 0, 0, Mkk, 0, 0);
    for (Index i = 0; i < n; ++i)
        Mkk(i, i) = Mkk(i, i) + shift;
}

void fill_offdiag_block(MatRealView Mij, Index ni, Index nj, Scalar scale)
{
    for (Index i = 0; i < ni; ++i)
        for (Index j = 0; j < nj; ++j)
            Mij(i, j) = scale * (::test::random() - 0.5);
}

// Build a known SPD block matrix in `M` using `sizes` and `edges`, with
// diagonal dominance guaranteed by a large shift on the diagonals.
void build_random_spd(BlockPdMatrix &M, const std::vector<Index> &sizes,
                      const std::vector<std::pair<Index, Index>> &edges)
{
    const Index N = static_cast<Index>(sizes.size());
    // Off-diagonal blocks first, with bounded magnitude.
    for (auto [a, b] : edges)
    {
        const Index i = std::max(a, b);
        const Index j = std::min(a, b);
        fill_offdiag_block(M.block(i, j), sizes[i], sizes[j], 0.1);
    }
    // Diagonals chosen large enough that the matrix is comfortably SPD even
    // for cliques of moderate size.
    for (Index k = 0; k < N; ++k)
    {
        const Scalar shift = static_cast<Scalar>(sizes[k]) + static_cast<Scalar>(N);
        fill_diagonal_block(M.block(k, k), sizes[k], shift);
    }
}
} // namespace

// --------------------------------------------------------------------------
// BlockSparsityPattern basic invariants.
// --------------------------------------------------------------------------
TEST(GraphSparsityTest, BuildAndQuery)
{
    std::vector<Index> sizes = {2, 3, 1, 4};
    std::vector<std::pair<Index, Index>> edges = {{0, 1}, {1, 2}, {0, 3}, {3, 0}, {2, 1}};
    BlockSparsityPattern sp(sizes, edges);

    EXPECT_EQ(sp.num_blocks(), 4);
    EXPECT_EQ(sp.total_size(), 2 + 3 + 1 + 4);
    EXPECT_EQ(sp.block_offset(0), 0);
    EXPECT_EQ(sp.block_offset(1), 2);
    EXPECT_EQ(sp.block_offset(4), 10);

    // Adjacency is deduped and sorted.
    EXPECT_EQ(sp.neighbors(0), (std::vector<Index>{1, 3}));
    EXPECT_EQ(sp.neighbors(1), (std::vector<Index>{0, 2}));
    EXPECT_EQ(sp.neighbors(2), (std::vector<Index>{1}));
    EXPECT_EQ(sp.neighbors(3), (std::vector<Index>{0}));

    // Column pattern in lower triangle: includes diagonal + rows > j.
    EXPECT_EQ(sp.column_pattern(0), (std::vector<Index>{0, 1, 3}));
    EXPECT_EQ(sp.column_pattern(1), (std::vector<Index>{1, 2}));
    EXPECT_EQ(sp.column_pattern(2), (std::vector<Index>{2}));
    EXPECT_EQ(sp.column_pattern(3), (std::vector<Index>{3}));
}

// --------------------------------------------------------------------------
// Elimination order is a valid permutation.
// --------------------------------------------------------------------------
TEST(GraphEliminationOrderTest, ValidPermutation)
{
    std::vector<Index> sizes(8, 2);
    std::vector<std::pair<Index, Index>> edges = {
        {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 7}, // chain
        {0, 7}                                                  // single long edge
    };
    BlockSparsityPattern sp(sizes, edges);
    BlockEliminationOrder ord(sp);
    std::vector<Index> seen(sp.num_blocks(), 0);
    for (Index k = 0; k < sp.num_blocks(); ++k)
    {
        EXPECT_GE(ord.perm()[k], 0);
        EXPECT_LT(ord.perm()[k], sp.num_blocks());
        seen[ord.perm()[k]]++;
    }
    for (Index k = 0; k < sp.num_blocks(); ++k)
    {
        EXPECT_EQ(seen[k], 1);
        EXPECT_EQ(ord.inv_perm()[ord.perm()[k]], k);
    }
}

// --------------------------------------------------------------------------
// Symbolic factorization: for a tridiagonal block chain the elimination tree
// is a path and there is no fill-in.
// --------------------------------------------------------------------------
TEST(GraphSymbolicTest, NoFillInOnChain)
{
    const Index N = 6;
    std::vector<Index> sizes(N, 3);
    std::vector<std::pair<Index, Index>> edges;
    for (Index k = 0; k < N - 1; ++k)
        edges.emplace_back(k, k + 1);
    BlockSparsityPattern sp(sizes, edges);
    BlockEliminationOrder ord(sp);
    BlockSymbolicFactorization sym(sp, ord);

    // Total non-diagonal entries in L equals the chain length, regardless of
    // ordering: a chain has no fill-in.
    Index total = 0;
    for (Index k = 0; k < N; ++k)
        total += static_cast<Index>(sym.lower_pattern(k).size());
    EXPECT_EQ(total, N - 1);
}

// --------------------------------------------------------------------------
// End-to-end correctness on a small chain (no fill-in expected).
// --------------------------------------------------------------------------
TEST(GraphCholeskySolverTest, SolveChain)
{
    const Index N = 5;
    std::vector<Index> sizes = {3, 2, 4, 2, 3};
    std::vector<std::pair<Index, Index>> edges;
    for (Index k = 0; k < N - 1; ++k)
        edges.emplace_back(k, k + 1);
    BlockSparsityPattern sp(sizes, edges);
    BlockPdMatrix M(sp);
    build_random_spd(M, sizes, edges);

    const Index m = sp.total_size();
    VecRealAllocated x_true(m);
    for (Index i = 0; i < m; ++i)
        x_true(i) = static_cast<Scalar>(i + 1) * 0.5;

    // b = -M * x_true so that M x = -b has the correct solution.
    VecRealAllocated zero(m);
    zero = 0.0;
    VecRealAllocated b(m);
    M.apply_on_right(x_true, 0.0, zero, b);
    vecsc(m, -1.0, b, 0);

    VecRealAllocated rhs_storage(m);
    VecRealView rhs_view(rhs_storage);
    veccp(m, b, 0, rhs_view, 0);

    LinearSystem<GraphType> ls(M, rhs_view);
    BlockCholeskySolver solver(sp);

    LinsolReturnFlag ret = solver.solve_in_place(ls);
    ASSERT_EQ(ret, LinsolReturnFlag::SUCCESS);

    // After solve_in_place, ls.rhs_view holds the solution x.
    for (Index i = 0; i < m; ++i)
        EXPECT_NEAR(rhs_view(i), x_true(i), 1e-8) << "row " << i;
}

// --------------------------------------------------------------------------
// Solve a non-trivial graph with fill-in: dense triangle on a subset of
// blocks. The minimum-degree ordering may or may not be optimal, but the
// numerical answer must match.
// --------------------------------------------------------------------------
TEST(GraphCholeskySolverTest, SolveWithFillIn)
{
    // 7 blocks; the first 4 form a clique, the last 3 form a chain hanging
    // off block 3.
    const Index N = 7;
    std::vector<Index> sizes = {2, 3, 2, 4, 3, 2, 2};
    std::vector<std::pair<Index, Index>> edges;
    for (Index i = 0; i < 4; ++i)
        for (Index j = i + 1; j < 4; ++j)
            edges.emplace_back(i, j);
    edges.emplace_back(3, 4);
    edges.emplace_back(4, 5);
    edges.emplace_back(5, 6);

    BlockSparsityPattern sp(sizes, edges);
    BlockPdMatrix M(sp);
    build_random_spd(M, sizes, edges);

    const Index m = sp.total_size();
    VecRealAllocated x_true(m);
    for (Index i = 0; i < m; ++i)
        x_true(i) = std::sin(0.7 * static_cast<Scalar>(i));

    VecRealAllocated zero(m);
    zero = 0.0;
    VecRealAllocated b(m);
    M.apply_on_right(x_true, 0.0, zero, b);
    vecsc(m, -1.0, b, 0);

    VecRealAllocated rhs_storage(m);
    VecRealView rhs_view(rhs_storage);
    veccp(m, b, 0, rhs_view, 0);

    LinearSystem<GraphType> ls(M, rhs_view);
    BlockCholeskySolver solver(sp);

    LinsolReturnFlag ret = solver.solve_in_place(ls);
    ASSERT_EQ(ret, LinsolReturnFlag::SUCCESS);

    for (Index i = 0; i < m; ++i)
        EXPECT_NEAR(rhs_view(i), x_true(i), 1e-7) << "row " << i;
}

// --------------------------------------------------------------------------
// solve_rhs should reuse a previously computed factorization.
// --------------------------------------------------------------------------
TEST(GraphCholeskySolverTest, SolveRhsReusesFactor)
{
    const Index N = 4;
    std::vector<Index> sizes = {2, 3, 2, 3};
    std::vector<std::pair<Index, Index>> edges = {{0, 1}, {1, 2}, {2, 3}, {0, 3}};
    BlockSparsityPattern sp(sizes, edges);
    BlockPdMatrix M(sp);
    build_random_spd(M, sizes, edges);

    const Index m = sp.total_size();
    VecRealAllocated x_true(m);
    for (Index i = 0; i < m; ++i)
        x_true(i) = 1.0 + 0.3 * static_cast<Scalar>(i);

    VecRealAllocated zero(m);
    zero = 0.0;
    VecRealAllocated b(m);
    M.apply_on_right(x_true, 0.0, zero, b);
    vecsc(m, -1.0, b, 0);

    VecRealAllocated rhs_storage(m);
    VecRealView rhs_view(rhs_storage);
    veccp(m, b, 0, rhs_view, 0);

    LinearSystem<GraphType> ls(M, rhs_view);
    BlockCholeskySolver solver(sp);

    ASSERT_EQ(solver.solve_in_place(ls), LinsolReturnFlag::SUCCESS);
    for (Index i = 0; i < m; ++i)
        EXPECT_NEAR(rhs_view(i), x_true(i), 1e-7);

    // Re-load a fresh RHS for a new x_true_2.
    VecRealAllocated x_true2(m);
    for (Index i = 0; i < m; ++i)
        x_true2(i) = std::cos(0.4 * static_cast<Scalar>(i));
    M.apply_on_right(x_true2, 0.0, zero, b);
    vecsc(m, -1.0, b, 0);
    veccp(m, b, 0, rhs_view, 0);

    ASSERT_EQ(solver.solve_in_place_rhs(ls), LinsolReturnFlag::SUCCESS);
    for (Index i = 0; i < m; ++i)
        EXPECT_NEAR(rhs_view(i), x_true2(i), 1e-7) << "row " << i;
}
