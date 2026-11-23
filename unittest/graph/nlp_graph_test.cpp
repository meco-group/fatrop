//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//
#include "fatrop/graph/graph_abstract.hpp"
#include "fatrop/graph/graph_blockwise_abstract.hpp"
#include "fatrop/graph/hessian.hpp"
#include "fatrop/graph/jacobian.hpp"
#include "fatrop/graph/nlp_graph.hpp"
#include "fatrop/graph/problem_info.hpp"
#include "fatrop/linear_algebra/blasfeo_operations.hpp"
#include "fatrop/linear_algebra/linear_algebra.hpp"

#include <gtest/gtest.h>
#include <memory>
#include <vector>

using namespace fatrop;

namespace
{
// Value written at entry (r, c) of a block tagged t.
Scalar tag_value(Index t, Index r, Index c) { return 100.0 * t + 10.0 * r + c; }

std::vector<Scalar> tagged_block(Index t, Index m, Index n)
{
    std::vector<Scalar> A(m * n);
    for (Index c = 0; c < n; ++c)
        for (Index r = 0; r < m; ++r)
            A[r + c * m] = tag_value(t, r, c);
    return A;
}

std::vector<Scalar> unpack(const MatRealView &A)
{
    std::vector<Scalar> out(A.m() * A.n());
    UNPACK_MAT(A.m(), A.n(), const_cast<MAT *>(&A.mat()), A.ai(), A.aj(), out.data(), A.m());
    return out;
}

// Mixed block sizes; edges deliberately not in fatrop's column-major storage order.
class TaggedGraph : public GraphAbstract
{
public:
    std::vector<Index> sizes{2, 3, 1, 2};
    std::vector<Index> ng{1, 0, 2, 0};
    std::vector<std::pair<Index, Index>> edges{{3, 0}, {2, 1}, {1, 0}, {3, 2}};

    GraphStructure get_structure() const override
    {
        GraphStructure st;
        st.block_sizes = sizes;
        st.ng_ineq = ng;
        st.off_diag_edges = edges;
        return st;
    }

    Index eval_hessian(const Scalar *, const Scalar *, const Scalar *,
                       GraphHessianView &H) override
    {
        for (Index k = 0; k < 4; ++k)
        {
            auto A = tagged_block(k, sizes[k], sizes[k]);
            MatRealView Hk = H.diag(k);
            pack_mat(sizes[k], sizes[k], A.data(), sizes[k], Hk, 0, 0);
        }
        for (Index e = 0; e < 4; ++e)
        {
            const Index m = sizes[edges[e].first], n = sizes[edges[e].second];
            auto A = tagged_block(10 + e, m, n);
            MatRealView He = H.off(e);
            pack_mat(m, n, A.data(), m, He, 0, 0);
        }
        return eval_gradient(nullptr, nullptr, H.grad());
    }

    Index eval_gradient(const Scalar *, const Scalar *, Scalar *grad) override
    {
        for (Index i = 0; i < 8; ++i)
            grad[i] += i;
        return 0;
    }

    Index eval_jacobian_ineq(const Scalar *, GraphJacobianView &J) override
    {
        for (Index k = 0; k < 4; ++k)
        {
            if (ng[k] == 0)
                continue;
            auto A = tagged_block(20 + k, sizes[k] + 1, ng[k]);
            MatRealView Jk = J.block(k);
            pack_mat(sizes[k] + 1, ng[k], A.data(), sizes[k] + 1, Jk, 0, 0);
        }
        return 0;
    }

    Index eval_gineq(const Scalar *, Scalar *g) override
    {
        for (Index i = 0; i < 3; ++i)
            g[i] = 1.0 + i;
        return 0;
    }

    Index eval_f(const Scalar *, const Scalar *, Scalar *res) override
    {
        *res = 0.0;
        return 0;
    }
    Index get_initial(Scalar *x) const override
    {
        for (Index i = 0; i < 8; ++i)
            x[i] = 0.0;
        return 0;
    }
    void apply_retraction(const Scalar *x, const Scalar *dx, const Scalar alpha,
                          Scalar *x_next) override
    {
        for (Index i = 0; i < 8; ++i)
            x_next[i] = x[i] + alpha * dx[i];
    }
};

// Same problem through the block-by-block interface; edges given in either
// orientation and with a duplicate.
class TaggedGraphBlockwise : public GraphBlockwiseAbstract
{
public:
    std::vector<Index> sizes{2, 3, 1, 2};
    std::vector<Index> ng{1, 0, 2, 0};

    Index get_num_blocks() const override { return 4; }
    Index get_block_size(Index k) const override { return sizes[k]; }
    Index get_ng_ineq(Index k) const override { return ng[k]; }
    std::vector<std::pair<Index, Index>> get_off_diag_edges() const override
    {
        return {{0, 3}, {2, 1}, {1, 0}, {3, 2}, {3, 0}};
    }

    Index eval_Hk(Index i, Index j, const Scalar *, const Scalar *, const Scalar *,
                  MAT *res) override
    {
        auto A = tagged_block(i == j ? i : 50 + 10 * i + j, sizes[i], sizes[j]);
        blasfeo_pack_mat_wrap(sizes[i], sizes[j], A.data(), sizes[i], res, 0, 0);
        return 0;
    }
    Index eval_grad_k(Index k, const Scalar *, const Scalar *, Scalar *res) override
    {
        Index off = 0;
        for (Index l = 0; l < k; ++l)
            off += sizes[l];
        for (Index a = 0; a < sizes[k]; ++a)
            res[a] = off + a;
        return 0;
    }
    Index eval_Ggt_ineq_k(Index k, const Scalar *, MAT *res) override
    {
        auto A = tagged_block(20 + k, sizes[k] + 1, ng[k]);
        blasfeo_pack_mat_wrap(sizes[k] + 1, ng[k], A.data(), sizes[k] + 1, res, 0, 0);
        return 0;
    }
    Index eval_gineq_k(Index k, const Scalar *, Scalar *res) override
    {
        for (Index a = 0; a < ng[k]; ++a)
            res[a] = (k == 0 ? 1.0 : 2.0) + a;
        return 0;
    }
    Index eval_f(const Scalar *, const Scalar *, Scalar *res) override
    {
        *res = 0.0;
        return 0;
    }
    Index get_initial(Scalar *x) const override
    {
        for (Index i = 0; i < 8; ++i)
            x[i] = 0.0;
        return 0;
    }
};

void expect_tagged(const MatRealView &A, Index t)
{
    const auto v = unpack(A);
    for (Index c = 0; c < A.n(); ++c)
        for (Index r = 0; r < A.m(); ++r)
            EXPECT_EQ(v[r + c * A.m()], tag_value(t, r, c)) << "tag " << t;
}
} // namespace

TEST(NlpGraphTest, HessianViewMapsUserEdges)
{
    auto graph = std::make_shared<TaggedGraph>();
    NlpGraph nlp(graph);
    const auto &dims = nlp.problem_dims();
    ProblemInfo<GraphProblem> info(dims);
    Hessian<GraphProblem> hess(dims);
    VecRealAllocated x(dims.nx), s(dims.ng_ineq_total), lam(dims.ng_ineq_total);
    x = 0.0;
    lam = 0.0;
    // Twice: the second call checks that fatrop zeroes before accumulating.
    for (int rep = 0; rep < 2; ++rep)
    {
        ASSERT_EQ(nlp.eval_lag_hess(info, 1.0, x, s, lam, hess), 0);
        for (Index k = 0; k < 4; ++k)
            expect_tagged(hess.H.block(k, k), k);
        for (Index e = 0; e < 4; ++e)
            expect_tagged(hess.H.block(graph->edges[e].first, graph->edges[e].second), 10 + e);
        for (Index i = 0; i < dims.nx_tan; ++i)
            EXPECT_EQ(hess.grad(i), i);
    }
}

TEST(NlpGraphTest, JacobianAndConstraints)
{
    auto graph = std::make_shared<TaggedGraph>();
    NlpGraph nlp(graph);
    const auto &dims = nlp.problem_dims();
    ProblemInfo<GraphProblem> info(dims);
    Jacobian<GraphProblem> jac(dims);
    VecRealAllocated x(dims.nx), s(dims.ng_ineq_total), res(info.number_of_eq_constraints);
    x = 0.0;
    s = 0.5;
    ASSERT_EQ(nlp.eval_constr_jac(info, x, s, jac), 0);
    expect_tagged(jac.Gg_ineqt[0], 20);
    expect_tagged(jac.Gg_ineqt[2], 22);
    ASSERT_EQ(nlp.eval_constraint_violation(info, x, s, res), 0);
    for (Index i = 0; i < 3; ++i)
        EXPECT_EQ(res(info.offset_g_eq_slack + i), 0.5 + i);
}

TEST(NlpGraphTest, RejectsUpperOrDuplicateEdges)
{
    auto upper = std::make_shared<TaggedGraph>();
    upper->edges[0] = {0, 3};
    EXPECT_THROW(NlpGraph{upper}, FatropException);
    auto dup = std::make_shared<TaggedGraph>();
    dup->edges[3] = {3, 0};
    EXPECT_THROW(NlpGraph{dup}, FatropException);
}

TEST(NlpGraphTest, BlockwiseAdapter)
{
    auto graph = std::make_shared<TaggedGraphBlockwise>();
    NlpGraph nlp(std::make_shared<GraphBlockwiseAdapter>(graph));
    const auto &dims = nlp.problem_dims();
    EXPECT_EQ(dims.off_diag_edges.size(), 4u);
    ProblemInfo<GraphProblem> info(dims);
    Hessian<GraphProblem> hess(dims);
    Jacobian<GraphProblem> jac(dims);
    VecRealAllocated x(dims.nx), s(dims.ng_ineq_total), lam(dims.ng_ineq_total),
        res(info.number_of_eq_constraints);
    x = 0.0;
    s = 0.5;
    lam = 0.0;
    ASSERT_EQ(nlp.eval_lag_hess(info, 1.0, x, s, lam, hess), 0);
    for (Index k = 0; k < 4; ++k)
        expect_tagged(hess.H.block(k, k), k);
    for (auto [i, j] : std::vector<std::pair<Index, Index>>{{3, 0}, {2, 1}, {1, 0}, {3, 2}})
        expect_tagged(hess.H.block(i, j), 50 + 10 * i + j);
    for (Index i = 0; i < dims.nx_tan; ++i)
        EXPECT_EQ(hess.grad(i), i);
    ASSERT_EQ(nlp.eval_constr_jac(info, x, s, jac), 0);
    expect_tagged(jac.Gg_ineqt[0], 20);
    expect_tagged(jac.Gg_ineqt[2], 22);
    ASSERT_EQ(nlp.eval_constraint_violation(info, x, s, res), 0);
    for (Index i = 0; i < 3; ++i)
        EXPECT_EQ(res(info.offset_g_eq_slack + i), 0.5 + i);
}
