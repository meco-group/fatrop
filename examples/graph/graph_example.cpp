//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//
// Simple graph-problem example: a quadratic objective coupled by a chain of
// neighbour interactions plus per-block box bounds. The problem is convex and
// has a unique solution that can be verified in closed form.
//
//   min   sum_k 0.5 * ||x_k - r_k||^2 + sum_{(i,j) in E} alpha * x_i^T x_j
//   s.t.  L_k <= x_k <= U_k                                      (block bounds)
//
// Each variable block has size 2; we use a chain of N = 4 blocks with edges
// (0,1), (1,2), (2,3). Bounds [-1, 1] are active for some components.
//
// This exercises:
//   - block-sparse Hessian assembly (diagonal + off-diagonal edges),
//   - block-local inequality constraints (box bounds via identity Jacobians),
//   - the BlockCholeskySolver in the augmented system,
//   - the full interior-point loop without restoration.
//

#include <fatrop/fatrop.hpp>
#include <fatrop/graph/graph_abstract.hpp>
#include <fatrop/graph/hessian.hpp>
#include <fatrop/graph/jacobian.hpp>
#include <fatrop/graph/nlp_graph.hpp>
#include <fatrop/graph/problem_info.hpp>
#include <fatrop/graph/problem_type.hpp>
#include <fatrop/ip_algorithm/ip_alg_builder.hpp>
#include <fatrop/ip_algorithm/ip_algorithm.hpp>
#include <fatrop/ip_algorithm/ip_data.hpp>

#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

using namespace fatrop;

namespace
{
constexpr Index N = 4;       // number of variable blocks
constexpr Index BLK = 2;     // each block has dimension 2
constexpr Scalar ALPHA = 0.3; // coupling strength
}

class GraphChainExample : public GraphAbstract
{
public:
    GraphChainExample()
    {
        r_.assign(N * BLK, 0.0);
        // Target points: r_k = (k + 1, -k).
        for (Index k = 0; k < N; ++k)
        {
            r_[k * BLK + 0] = static_cast<Scalar>(k + 1);
            r_[k * BLK + 1] = -static_cast<Scalar>(k);
        }
    }

    Index get_num_blocks() const override { return N; }
    Index get_block_size(Index) const override { return BLK; }
    Index get_ng_ineq(Index) const override { return BLK; }

    std::vector<std::pair<Index, Index>> get_off_diag_edges() const override
    {
        std::vector<std::pair<Index, Index>> edges;
        for (Index k = 0; k + 1 < N; ++k)
            edges.emplace_back(k, k + 1);
        return edges;
    }

    Index eval_Hk(Index i, Index j, const Scalar *objective_scale, const Scalar * /*x*/,
                  const Scalar * /*mult*/, MAT *res) override
    {
        const Scalar s = objective_scale[0];
        blasfeo_gese_wrap(res->m, res->n, 0.0, res, 0, 0);
        if (i == j)
        {
            // Hessian diagonal block: identity (from 0.5 * ||x_k - r_k||^2).
            for (Index a = 0; a < BLK; ++a)
                blasfeo_matel_wrap(res, a, a) = s * 1.0;
        }
        else
        {
            // Off-diagonal block (i > j): alpha * I  (from alpha * x_i^T x_j).
            for (Index a = 0; a < BLK; ++a)
                blasfeo_matel_wrap(res, a, a) = s * ALPHA;
        }
        return 0;
    }

    Index eval_grad_k(Index k, const Scalar *objective_scale, const Scalar *x, Scalar *res) override
    {
        const Scalar s = objective_scale[0];
        // grad_k = (x_k - r_k) + alpha * (sum of neighbour x's).
        for (Index a = 0; a < BLK; ++a)
            res[a] = s * (x[k * BLK + a] - r_[k * BLK + a]);
        // Neighbours: (k-1) and (k+1) in the chain.
        if (k > 0)
            for (Index a = 0; a < BLK; ++a)
                res[a] += s * ALPHA * x[(k - 1) * BLK + a];
        if (k + 1 < N)
            for (Index a = 0; a < BLK; ++a)
                res[a] += s * ALPHA * x[(k + 1) * BLK + a];
        return 0;
    }

    Index eval_Ggt_ineq_k(Index k, const Scalar *x, MAT *res) override
    {
        // Inequality g_ineq_k(x_k) = x_k, so the (BLK+1) x BLK transposed
        // Jacobian's top part is the identity. The trailing row stores the
        // current value of g_ineq_k.
        blasfeo_gese_wrap(res->m, res->n, 0.0, res, 0, 0);
        for (Index a = 0; a < BLK; ++a)
        {
            blasfeo_matel_wrap(res, a, a) = 1.0;
            blasfeo_matel_wrap(res, BLK, a) = x[k * BLK + a];
        }
        return 0;
    }

    Index eval_gineq_k(Index k, const Scalar *x, Scalar *res) override
    {
        for (Index a = 0; a < BLK; ++a)
            res[a] = x[k * BLK + a];
        return 0;
    }

    Index eval_f(const Scalar *objective_scale, const Scalar *x, Scalar *res) override
    {
        const Scalar s = objective_scale[0];
        Scalar v = 0.0;
        for (Index k = 0; k < N; ++k)
        {
            for (Index a = 0; a < BLK; ++a)
            {
                Scalar d = x[k * BLK + a] - r_[k * BLK + a];
                v += 0.5 * d * d;
            }
        }
        for (Index k = 0; k + 1 < N; ++k)
            for (Index a = 0; a < BLK; ++a)
                v += ALPHA * x[k * BLK + a] * x[(k + 1) * BLK + a];
        *res = s * v;
        return 0;
    }

    Index get_bounds(Scalar *lower, Scalar *upper) const override
    {
        for (Index k = 0; k < N; ++k)
            for (Index a = 0; a < BLK; ++a)
            {
                lower[k * BLK + a] = -1.0;
                upper[k * BLK + a] = 1.0;
            }
        return 0;
    }

    Index get_initial(Scalar *x) const override
    {
        for (Index i = 0; i < N * BLK; ++i)
            x[i] = 0.0;
        return 0;
    }

private:
    std::vector<Scalar> r_;
};

int main()
{
    auto problem = std::make_shared<GraphChainExample>();
    OptionRegistry opts;
    IpAlgBuilder<GraphProblem> builder(std::make_shared<NlpGraph>(problem));
    auto alg = builder.with_options_registry(&opts).build();

    Timer t;
    t.start();
    IpSolverReturnFlag ret = alg->optimize();
    const auto elapsed = t.stop();
    auto data = builder.get_ipdata();

    std::cout << "\n=== fatrop graph-problem IP ===\n";
    std::cout << "Elapsed: " << elapsed << "\n";
    std::cout << data->timing_statistics() << "\n";
    std::cout << "return flag = " << int(ret) << "\n";
    const VecRealView &x = data->current_iterate().primal_x();
    std::cout << "x* = (";
    for (Index i = 0; i < x.m(); ++i)
        std::cout << x(i) << (i + 1 == x.m() ? "" : ", ");
    std::cout << ")\n";

    return ret == IpSolverReturnFlag::Success ? 0 : 1;
}
