//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//
#include "fatrop/graph/block_sparsity.hpp"
#include "fatrop/graph/elimination_order.hpp"
#include "fatrop/graph/symbolic_factorization.hpp"

#include <algorithm>
#include <gtest/gtest.h>
#include <numeric>
#include <random>
#include <set>
#include <vector>

using namespace fatrop;

namespace
{
using Edges = std::vector<std::pair<Index, Index>>;

void expect_valid(const BlockSparsityPattern &sp, const BlockEliminationOrder &ord)
{
    const Index n = sp.num_blocks();
    ASSERT_EQ(static_cast<Index>(ord.perm().size()), n);
    ASSERT_EQ(static_cast<Index>(ord.inv_perm().size()), n);
    std::vector<Index> seen(n, 0);
    for (Index k = 0; k < n; ++k)
    {
        ASSERT_GE(ord.perm()[k], 0);
        ASSERT_LT(ord.perm()[k], n);
        seen[ord.perm()[k]]++;
        EXPECT_EQ(ord.inv_perm()[ord.perm()[k]], k);
    }
    for (Index k = 0; k < n; ++k)
        EXPECT_EQ(seen[k], 1);
}

// Scalar nonzeros of the strictly-lower block part of L for elimination order
// perm, by explicit elimination.
Index lower_nnz(const BlockSparsityPattern &sp, const std::vector<Index> &perm)
{
    const Index n = sp.num_blocks();
    std::vector<std::set<Index>> adj(n);
    for (Index v = 0; v < n; ++v)
        adj[v].insert(sp.neighbors(v).begin(), sp.neighbors(v).end());
    Index nnz = 0;
    for (Index v : perm)
    {
        std::vector<Index> nb(adj[v].begin(), adj[v].end());
        for (Index u : nb)
        {
            nnz += sp.block_size(u) * sp.block_size(v);
            adj[u].erase(v);
            for (Index t : nb)
                if (t != u)
                    adj[u].insert(t);
        }
    }
    return nnz;
}

Index lower_nnz(const BlockSparsityPattern &sp, const BlockSymbolicFactorization &sym,
                const BlockEliminationOrder &ord)
{
    Index nnz = 0;
    for (Index j = 0; j < sp.num_blocks(); ++j)
        for (Index i : sym.lower_pattern(j))
            nnz += sp.block_size(ord.perm()[i]) * sp.block_size(ord.perm()[j]);
    return nnz;
}

Index natural_nnz(const BlockSparsityPattern &sp)
{
    std::vector<Index> perm(sp.num_blocks());
    std::iota(perm.begin(), perm.end(), 0);
    return lower_nnz(sp, perm);
}

Edges grid(Index nx, Index ny)
{
    Edges e;
    for (Index x = 0; x < nx; ++x)
        for (Index y = 0; y < ny; ++y)
        {
            if (x + 1 < nx)
                e.emplace_back(x * ny + y, (x + 1) * ny + y);
            if (y + 1 < ny)
                e.emplace_back(x * ny + y, x * ny + y + 1);
        }
    return e;
}

Edges random_graph(Index n, Index m, unsigned seed)
{
    std::mt19937 rng(seed);
    std::uniform_int_distribution<Index> u(0, n - 1);
    Edges e;
    for (Index k = 1; k < n; ++k) // spanning chain plus random loop closures
        e.emplace_back(k - 1, k);
    for (Index k = 0; k < m; ++k)
    {
        const Index a = u(rng), b = u(rng);
        if (a != b)
            e.emplace_back(a, b);
    }
    return e;
}

std::vector<Index> random_sizes(Index n, Index lo, Index hi, unsigned seed)
{
    std::mt19937 rng(seed);
    std::uniform_int_distribution<Index> u(lo, hi);
    std::vector<Index> s(n);
    for (auto &x : s)
        x = u(rng);
    return s;
}

void check(const std::vector<Index> &sizes, const Edges &edges, bool no_worse_than_natural)
{
    BlockSparsityPattern sp(sizes, edges);
    BlockEliminationOrder ord(sp);
    expect_valid(sp, ord);
    BlockEliminationOrder again(sp);
    EXPECT_EQ(ord.perm(), again.perm()); // deterministic
    BlockSymbolicFactorization sym(sp, ord);
    const Index nnz = lower_nnz(sp, ord.perm());
    EXPECT_EQ(lower_nnz(sp, sym, ord), nnz);
    if (no_worse_than_natural)
        EXPECT_LE(nnz, natural_nnz(sp));
}
} // namespace

TEST(GraphAmdTest, SingleBlock) { check({4}, {}, true); }

TEST(GraphAmdTest, NoEdges) { check({1, 2, 3, 0, 5}, {}, true); }

TEST(GraphAmdTest, Disconnected)
{
    // Two cycles, a path, and isolated blocks.
    check({2, 3, 1, 2, 4, 2, 3, 1, 1, 6, 2},
          {{0, 1}, {1, 2}, {2, 0}, {3, 4}, {4, 5}, {5, 6}, {6, 3}, {7, 8}}, true);
}

TEST(GraphAmdTest, DenseClique)
{
    const Index n = 12;
    Edges e;
    for (Index i = 0; i < n; ++i)
        for (Index j = 0; j < i; ++j)
            e.emplace_back(i, j);
    check(random_sizes(n, 1, 6, 1), e, true);
    // All blocks of a clique are indistinguishable: one supernode.
    BlockSparsityPattern sp(std::vector<Index>(n, 3), e);
    BlockEliminationOrder ord(sp);
    BlockSymbolicFactorization sym(sp, ord);
    EXPECT_EQ(sym.num_supernodes(), 1);
}

TEST(GraphAmdTest, ZeroSizeBlocks)
{
    std::vector<Index> sizes = random_sizes(60, 0, 4, 2);
    sizes[0] = 0;
    check(sizes, random_graph(60, 40, 3), false);
    check(std::vector<Index>(10, 0), random_graph(10, 8, 4), false);
}

TEST(GraphAmdTest, Arrow)
{
    // Hub first: the natural order fills in completely.
    const Index n = 30;
    Edges e;
    for (Index k = 1; k < n; ++k)
        e.emplace_back(0, k);
    BlockSparsityPattern sp(std::vector<Index>(n, 2), e);
    BlockEliminationOrder ord(sp);
    expect_valid(sp, ord);
    EXPECT_EQ(lower_nnz(sp, ord.perm()), (n - 1) * 4);
}

TEST(GraphAmdTest, Grid) { check(random_sizes(15 * 15, 1, 6, 5), grid(15, 15), true); }

TEST(GraphAmdTest, PoseGraphLike)
{
    for (unsigned seed = 0; seed < 5; ++seed)
        check(random_sizes(400, 3, 6, 10 + seed), random_graph(400, 120, 20 + seed), true);
}
