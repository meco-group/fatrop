#ifndef FATROP_EXAMPLES_LIMITED_MEMORY_PROBLEM_HPP
#define FATROP_EXAMPLES_LIMITED_MEMORY_PROBLEM_HPP

#include <algorithm>
#include <cmath>
#include <fatrop/fatrop.hpp>
#include <fatrop/graph/graph_abstract.hpp>
#include <fatrop/graph/hessian.hpp>
#include <fatrop/graph/jacobian.hpp>
#include <fatrop/graph/nlp_graph.hpp>
#include <fatrop/graph/problem_info.hpp>
#include <fatrop/graph/problem_type.hpp>
#include <stdexcept>
#include <vector>

namespace limited_memory_example
{
using namespace fatrop;

// Shared by Fatrop and Ipopt. With equality=true:
//   c(x) = x[0]^2 + sum_{i>0} x[i] - n = 0,
//   f(x) = sum_i ((i+1)/2 * (x[i]-1)^2 + (x[i]-1)^4/4) - c(x).
// On c=0 the objective is nonnegative and x*=1 is its unique minimizer;
// the equality multiplier is 1, so the Jacobian term in BFGS is essential.
// Without the equality, add chain coupling instead of -c; x*=1 still holds.
// Both cases have nonlinear inequalities 0 <= x[i] + 0.1*x[i]^2 <= 2.4.
struct Model
{
    static constexpr Index n = 8;
    bool equality = true;
    Scalar constraint(const Scalar *x) const
    {
        Scalar c = x[0] * x[0] - n;
        for (Index i = 1; i < n; ++i)
            c += x[i];
        return c;
    }
    Scalar objective(const Scalar *x) const
    {
        Scalar f = equality ? -constraint(x) : 0;
        for (Index i = 0; i < n; ++i)
        {
            const Scalar e = x[i] - 1;
            f += 0.5 * (i + 1) * e * e + 0.25 * e * e * e * e;
            if (!equality && i + 1 < n)
                f += 0.2 * (x[i] - x[i + 1]) * (x[i] - x[i + 1]);
        }
        return f;
    }
    Scalar gradient(Index i, const Scalar *x) const
    {
        const Scalar e = x[i] - 1;
        Scalar g = (i + 1) * e + e * e * e;
        if (equality)
            g -= i == 0 ? 2 * x[0] : 1;
        else
        {
            if (i > 0) g += 0.4 * (x[i] - x[i - 1]);
            if (i + 1 < n) g += 0.4 * (x[i] - x[i + 1]);
        }
        return g;
    }
    static Scalar inequality(Scalar x) { return x + 0.1 * x * x; }
    static Scalar inequality_derivative(Scalar x) { return 1 + 0.2 * x; }
    static void initial(Scalar *x)
    {
        for (Index i = 0; i < n; ++i) x[i] = 0.3 + 0.12 * i;
    }
};

class DenseProblem : public DenseAbstract
{
public:
    explicit DenseProblem(bool equality = true) : model{equality} {}
    Model model;
    Index hessian_calls = 0;
    Index get_nx() const override { return Model::n; }
    Index get_ng() const override { return model.equality ? 1 : 0; }
    Index get_ng_ineq() const override { return Model::n; }
    Index eval_Hh(const Scalar *, const Scalar *, const Scalar *, MAT *) override
    {
        ++hessian_calls;
        throw std::runtime_error("L-BFGS must not evaluate the exact Hessian");
    }
    Index eval_Ggt(const Scalar *x, MAT *res) override
    {
        if (!model.equality) return 0;
        for (Index i = 0; i < Model::n; ++i)
            blasfeo_matel_wrap(res, i, 0) = i == 0 ? 2 * x[0] : 1;
        blasfeo_matel_wrap(res, Model::n, 0) = model.constraint(x);
        return 0;
    }
    Index eval_Ggt_ineq(const Scalar *x, MAT *res) override
    {
        blasfeo_gese_wrap(res->m, res->n, 0, res, 0, 0);
        for (Index i = 0; i < Model::n; ++i)
        {
            blasfeo_matel_wrap(res, i, i) = Model::inequality_derivative(x[i]);
            blasfeo_matel_wrap(res, Model::n, i) = Model::inequality(x[i]);
        }
        return 0;
    }
    Index eval_g(const Scalar *x, Scalar *res) override
    {
        if (model.equality) res[0] = model.constraint(x);
        return 0;
    }
    Index eval_gineq(const Scalar *x, Scalar *res) override
    {
        for (Index i = 0; i < Model::n; ++i) res[i] = Model::inequality(x[i]);
        return 0;
    }
    Index eval_grad(const Scalar *scale, const Scalar *x, Scalar *res) override
    {
        for (Index i = 0; i < Model::n; ++i) res[i] = *scale * model.gradient(i, x);
        return 0;
    }
    Index eval_f(const Scalar *scale, const Scalar *x, Scalar *res) override
    {
        *res = *scale * model.objective(x);
        return 0;
    }
    Index get_bounds(Scalar *lower, Scalar *upper) const override
    {
        for (Index i = 0; i < Model::n; ++i) { lower[i] = 0; upper[i] = 2.4; }
        return 0;
    }
    Index get_initial(Scalar *x) const override { Model::initial(x); return 0; }
};

// A one-stage OCP exposes exactly the dense problem, including multiplier order.
class OcpProblem : public OcpAbstract
{
public:
    DenseProblem dense;
    Index get_nx(Index) const override { return Model::n; }
    Index get_nu(Index) const override { return 0; }
    Index get_ng(Index) const override { return 1; }
    Index get_ng_ineq(Index) const override { return Model::n; }
    Index get_horizon_length() const override { return 1; }
    Index eval_BAbt(const Scalar *, const Scalar *, const Scalar *, MAT *, Index) override { return 0; }
    Index eval_b(const Scalar *, const Scalar *, const Scalar *, Scalar *, Index) override { return 0; }
    Index eval_RSQrqt(const Scalar *a, const Scalar *, const Scalar *x, const Scalar *,
                     const Scalar *lam, const Scalar *, MAT *r, Index) override
    { return dense.eval_Hh(a, x, lam, r); }
    Index eval_Ggt(const Scalar *, const Scalar *x, MAT *r, Index) override
    { return dense.eval_Ggt(x, r); }
    Index eval_Ggt_ineq(const Scalar *, const Scalar *x, MAT *r, Index) override
    { return dense.eval_Ggt_ineq(x, r); }
    Index eval_g(const Scalar *, const Scalar *x, Scalar *r, Index) override
    { return dense.eval_g(x, r); }
    Index eval_gineq(const Scalar *, const Scalar *x, Scalar *r, Index) override
    { return dense.eval_gineq(x, r); }
    Index eval_rq(const Scalar *a, const Scalar *, const Scalar *x, Scalar *r, Index) override
    { return dense.eval_grad(a, x, r); }
    Index eval_L(const Scalar *a, const Scalar *, const Scalar *x, Scalar *r, Index) override
    { return dense.eval_f(a, x, r); }
    Index get_bounds(Scalar *l, Scalar *u, Index) const override { return dense.get_bounds(l, u); }
    Index get_initial_xk(Scalar *x, Index) const override { return dense.get_initial(x); }
    Index get_initial_uk(Scalar *, Index) const override { return 0; }
};

class GraphProblem : public GraphAbstract
{
public:
    Model model{false};
    Index hessian_calls = 0;
    Index get_num_blocks() const override { return Model::n; }
    Index get_block_size(Index) const override { return 1; }
    Index get_ng_ineq(Index) const override { return 1; }
    std::vector<std::pair<Index, Index>> get_off_diag_edges() const override
    {
        std::vector<std::pair<Index, Index>> edges;
        for (Index i = 1; i < Model::n; ++i) edges.emplace_back(i, i - 1);
        return edges;
    }
    Index eval_Hk(Index, Index, const Scalar *, const Scalar *, const Scalar *, MAT *) override
    {
        ++hessian_calls;
        throw std::runtime_error("L-BFGS must not evaluate graph Hessian blocks");
    }
    Index eval_grad_k(Index k, const Scalar *a, const Scalar *x, Scalar *r) override
    { *r = *a * model.gradient(k, x); return 0; }
    Index eval_Ggt_ineq_k(Index k, const Scalar *x, MAT *r) override
    {
        blasfeo_matel_wrap(r, 0, 0) = Model::inequality_derivative(x[k]);
        blasfeo_matel_wrap(r, 1, 0) = Model::inequality(x[k]);
        return 0;
    }
    Index eval_gineq_k(Index k, const Scalar *x, Scalar *r) override
    { *r = Model::inequality(x[k]); return 0; }
    Index eval_f(const Scalar *a, const Scalar *x, Scalar *r) override
    { *r = *a * model.objective(x); return 0; }
    Index get_bounds(Scalar *l, Scalar *u) const override
    { for (Index i = 0; i < Model::n; ++i) { l[i] = 0; u[i] = 2.4; } return 0; }
    Index get_initial(Scalar *x) const override { Model::initial(x); return 0; }
};

// Six stages with nonlinear dynamics, a nonlinear initial equality, and bounds.
// Subtracting the sum of the equalities from a positive objective makes every
// dynamics multiplier and the initial-condition multiplier nonzero at x*=1,u*=0.
class DynamicsProblem : public OcpAbstract
{
public:
    static constexpr Index horizon = 6;
    Index hessian_calls = 0;
    Index get_nx(Index) const override { return 1; }
    Index get_nu(Index k) const override { return k + 1 < horizon ? 1 : 0; }
    Index get_ng(Index k) const override { return k == 0 ? 1 : 0; }
    Index get_ng_ineq(Index k) const override { return 1 + get_nu(k); }
    Index get_horizon_length() const override { return horizon; }
    Index eval_BAbt(const Scalar *xn, const Scalar *u, const Scalar *x, MAT *r, Index) override
    {
        blasfeo_matel_wrap(r, 0, 0) = 1 + 0.2 * u[0];
        blasfeo_matel_wrap(r, 1, 0) = 1;
        blasfeo_matel_wrap(r, 2, 0) = x[0] + u[0] + 0.1 * u[0] * u[0] - xn[0];
        return 0;
    }
    Index eval_b(const Scalar *xn, const Scalar *u, const Scalar *x, Scalar *r, Index) override
    { r[0] = x[0] + u[0] + 0.1 * u[0] * u[0] - xn[0]; return 0; }
    Index eval_RSQrqt(const Scalar *, const Scalar *, const Scalar *, const Scalar *,
                     const Scalar *, const Scalar *, MAT *, Index) override
    { ++hessian_calls; throw std::runtime_error("Multistage L-BFGS requested a Hessian"); }
    Index eval_Ggt(const Scalar *, const Scalar *x, MAT *r, Index k) override
    {
        if (k == 0)
        {
            blasfeo_matel_wrap(r, 0, 0) = 0;
            blasfeo_matel_wrap(r, 1, 0) = 2 * x[0];
            blasfeo_matel_wrap(r, 2, 0) = x[0] * x[0] - 1;
        }
        return 0;
    }
    Index eval_g(const Scalar *, const Scalar *x, Scalar *r, Index k) override
    { if (k == 0) r[0] = x[0] * x[0] - 1; return 0; }
    Index eval_Ggt_ineq(const Scalar *u, const Scalar *x, MAT *r, Index k) override
    {
        const Index nu = get_nu(k);
        blasfeo_gese_wrap(r->m, r->n, 0, r, 0, 0);
        if (nu)
        {
            blasfeo_matel_wrap(r, 0, 0) = 1;
            blasfeo_matel_wrap(r, 2, 0) = u[0];
        }
        blasfeo_matel_wrap(r, nu, nu) = Model::inequality_derivative(x[0]);
        blasfeo_matel_wrap(r, nu + 1, nu) = Model::inequality(x[0]);
        return 0;
    }
    Index eval_gineq(const Scalar *u, const Scalar *x, Scalar *r, Index k) override
    {
        const Index nu = get_nu(k);
        if (nu) r[0] = u[0];
        r[nu] = Model::inequality(x[0]);
        return 0;
    }
    Index eval_rq(const Scalar *a, const Scalar *u, const Scalar *x, Scalar *r, Index k) override
    {
        const Index nu = get_nu(k);
        if (nu) r[0] = *a * (u[0] + u[0] * u[0] * u[0] - 1 - 0.2 * u[0]);
        const Scalar e = x[0] - 1;
        Scalar g = (k + 1) * e + e * e * e;
        if (k == 0) g -= 1 + 2 * x[0];
        if (k + 1 == horizon) g += 1;
        r[nu] = *a * g;
        return 0;
    }
    Index eval_L(const Scalar *a, const Scalar *u, const Scalar *x, Scalar *r, Index k) override
    {
        const Scalar e = x[0] - 1;
        Scalar f = 0.5 * (k + 1) * e * e + 0.25 * e * e * e * e;
        if (get_nu(k))
            f += 0.5 * u[0] * u[0] + 0.25 * std::pow(u[0], 4) - u[0] - 0.1 * u[0] * u[0];
        if (k == 0) f -= x[0] + x[0] * x[0] - 1;
        if (k + 1 == horizon) f += x[0];
        *r = *a * f;
        return 0;
    }
    Index get_bounds(Scalar *l, Scalar *u, Index k) const override
    {
        const Index nu = get_nu(k);
        if (nu) { l[0] = -1; u[0] = 1; }
        l[nu] = 0; u[nu] = 2.4;
        return 0;
    }
    Index get_initial_xk(Scalar *x, Index) const override { x[0] = 0.8; return 0; }
    Index get_initial_uk(Scalar *u, Index k) const override
    { if (get_nu(k)) u[0] = 0.2; return 0; }
};

struct Result
{
    bool success = false;
    Index iterations = 0;
    Index hessian_calls = 0;
    Scalar objective = 0;
    Scalar max_error = 0;
    Scalar feasibility = 0;
    std::vector<Scalar> x;
};

inline void evaluate_result(Result &r, const Model &model)
{
    r.objective = model.objective(r.x.data());
    r.feasibility = model.equality ? std::abs(model.constraint(r.x.data())) : 0;
    for (Scalar x : r.x)
    {
        r.max_error = std::max(r.max_error, std::abs(x - 1));
        const Scalar d = Model::inequality(x);
        r.feasibility = std::max(r.feasibility, std::max(-d, d - 2.4));
    }
}

template <typename ProblemType>
Result solve_fatrop(std::shared_ptr<Nlp<ProblemType>> nlp, const Model &model,
                    Index history = 6, const std::string &approximation = "limited-memory")
{
    OptionRegistry options;
    IpAlgBuilder<ProblemType> builder(nlp);
    auto algorithm = builder.with_options_registry(&options).build();
    options.set_option("hessian_approximation", approximation);
    options.set_option("limited_memory_max_history", history);
    options.set_option("tolerance", Scalar(1e-8));
    options.set_option("constr_viol_tol", Scalar(1e-8));
    options.set_option("compl_inf_tol", Scalar(1e-8));
    options.set_option("tol_acceptable", Scalar(1e-8));
    options.set_option("max_iter", Index(300));
    options.set_option("print_level", Index(0));
    Result result;
    result.success = algorithm->optimize() == IpSolverReturnFlag::Success;
    auto data = builder.get_ipdata();
    result.iterations = data->iteration_number();
    const auto &x = data->current_iterate().primal_x();
    for (Index i = 0; i < x.m(); ++i) result.x.push_back(x(i));
    evaluate_result(result, model);
    return result;
}
} // namespace limited_memory_example
#endif
