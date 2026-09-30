#include <fatrop/fatrop.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <vector>

using namespace fatrop;

namespace
{
// Coupled OCP with configurable state dimension.
class LargeStateOcp : public OcpAbstract
{
public:
    explicit LargeStateOcp(Index states, Index stages = 60)
        : states_(states), stages_(stages) {}
    Index get_nx(Index) const override { return states_; }
    Index get_nu(Index k) const override { return k + 1 < stages_ ? 8 : 0; }
    Index get_ng(Index k) const override { return k == 0 ? states_ : 0; }
    Index get_ng_ineq(Index k) const override { return get_nu(k); }
    Index get_horizon_length() const override { return stages_; }

    Index eval_BAbt(const Scalar *next, const Scalar *u, const Scalar *x,
                    MAT *r, Index) override
    {
        blasfeo_gese_wrap(r->m, r->n, 0., r, 0, 0);
        for (Index i = 0; i < states_; ++i)
        {
            const Index j = i % 8;
            blasfeo_matel_wrap(r, j, i) = 0.15 + 0.02 * u[j];
            blasfeo_matel_wrap(r, 8 + i, i) = 0.9;
            blasfeo_matel_wrap(r, 8 + states_, i) =
                0.9 * x[i] + 0.15 * u[j] + 0.01 * u[j] * u[j] - next[i];
        }
        return 0;
    }
    Index eval_b(const Scalar *next, const Scalar *u, const Scalar *x,
                 Scalar *r, Index) override
    {
        for (Index i = 0; i < states_; ++i)
        {
            const Scalar v = u[i % 8];
            r[i] = 0.9 * x[i] + 0.15 * v + 0.01 * v * v - next[i];
        }
        return 0;
    }
    Index eval_RSQrqt(const Scalar *scale, const Scalar *, const Scalar *x,
                      const Scalar *lam_dyn, const Scalar *, const Scalar *,
                      MAT *r, Index k) override
    {
        blasfeo_gese_wrap(r->m, r->n, 0., r, 0, 0);
        const Index nu = get_nu(k);
        for (Index j = 0; j < nu; ++j)
        {
            Scalar curvature = 0.2 * *scale;
            for (Index i = j; i < states_; i += 8)
                curvature += 0.02 * lam_dyn[i];
            blasfeo_matel_wrap(r, j, j) = curvature;
        }
        for (Index i = 0; i < states_; ++i)
            blasfeo_matel_wrap(r, nu + i, nu + i) = 2. * *scale;
        return 0;
    }
    Index eval_Ggt(const Scalar *, const Scalar *x, MAT *r, Index k) override
    {
        if (k != 0) return 0;
        blasfeo_gese_wrap(r->m, r->n, 0., r, 0, 0);
        for (Index i = 0; i < states_; ++i)
        {
            blasfeo_matel_wrap(r, 8 + i, i) = 1.;
            blasfeo_matel_wrap(r, 8 + states_, i) = x[i];
        }
        return 0;
    }
    Index eval_g(const Scalar *, const Scalar *x, Scalar *r, Index k) override
    {
        if (k == 0) for (Index i = 0; i < states_; ++i) r[i] = x[i];
        return 0;
    }
    Index eval_Ggt_ineq(const Scalar *u, const Scalar *, MAT *r, Index k) override
    {
        if (get_nu(k) == 0) return 0;
        blasfeo_gese_wrap(r->m, r->n, 0., r, 0, 0);
        for (Index j = 0; j < 8; ++j)
        {
            blasfeo_matel_wrap(r, j, j) = 1.;
            blasfeo_matel_wrap(r, 8 + states_, j) = u[j];
        }
        return 0;
    }
    Index eval_gineq(const Scalar *u, const Scalar *, Scalar *r, Index k) override
    {
        if (get_nu(k)) for (Index j = 0; j < 8; ++j) r[j] = u[j];
        return 0;
    }
    Index eval_rq(const Scalar *scale, const Scalar *u, const Scalar *x,
                  Scalar *r, Index k) override
    {
        const Index nu = get_nu(k);
        for (Index j = 0; j < nu; ++j) r[j] = 0.2 * *scale * u[j];
        for (Index i = 0; i < states_; ++i)
            r[nu + i] = 2. * *scale * (x[i] - 1.);
        return 0;
    }
    Index eval_L(const Scalar *scale, const Scalar *u, const Scalar *x,
                 Scalar *r, Index k) override
    {
        Scalar f = 0.;
        for (Index j = 0; j < get_nu(k); ++j) f += 0.1 * u[j] * u[j];
        for (Index i = 0; i < states_; ++i)
        {
            const Scalar e = x[i] - 1.;
            f += e * e;
        }
        *r = *scale * f;
        return 0;
    }
    Index get_bounds(Scalar *lower, Scalar *upper, Index k) const override
    {
        for (Index j = 0; j < get_nu(k); ++j)
        {
            lower[j] = -2.;
            upper[j] = 2.;
        }
        return 0;
    }
    Index get_initial_xk(Scalar *x, Index) const override
    { for (Index i = 0; i < states_; ++i) x[i] = 0.; return 0; }
    Index get_initial_uk(Scalar *u, Index k) const override
    { for (Index j = 0; j < get_nu(k); ++j) u[j] = 0.; return 0; }

private:
    Index states_, stages_;
};

struct Runner
{
    OptionRegistry options;
    IpAlgBuilder<OcpType> builder;
    std::shared_ptr<IpAlgorithm<OcpType>> algorithm;
    std::vector<double> time_per_iteration_us;
    std::vector<double> total_us;
    Index iterations = 0;

    Runner(Index states, const char *mode)
        : builder(std::make_shared<NlpOcp>(std::make_shared<LargeStateOcp>(states)))
    {
        algorithm = builder.with_options_registry(&options).build();
        options.set_option("hessian_approximation", std::string(mode));
        options.set_option("tolerance", Scalar(1e-8));
        options.set_option("constr_viol_tol", Scalar(1e-8));
        options.set_option("compl_inf_tol", Scalar(1e-8));
        options.set_option("tol_acceptable", Scalar(1e-8));
        options.set_option("max_iter", Index(300));
        options.set_option("print_level", Index(0));
        options.set_option("suppress_banner", true);
    }
    void run(bool record)
    {
        const auto start = std::chrono::steady_clock::now();
        const auto status = algorithm->optimize();
        const double elapsed = std::chrono::duration<double, std::micro>(
            std::chrono::steady_clock::now() - start).count();
        if (status != IpSolverReturnFlag::Success)
        {
            std::fprintf(stderr, "solve failed: %d\n", int(status));
            std::abort();
        }
        iterations = builder.get_ipdata()->iteration_number();
        if (record)
        {
            total_us.push_back(elapsed);
            time_per_iteration_us.push_back(elapsed / iterations);
        }
    }
    void print(Index states, const char *mode)
    {
        std::sort(total_us.begin(), total_us.end());
        std::sort(time_per_iteration_us.begin(), time_per_iteration_us.end());
        const auto middle = total_us.size() / 2;
        const auto data = builder.get_ipdata();
        const auto &info = data->info();
        const auto &solution = data->current_iterate().primal_x();
        const auto &timings = data->timing_statistics();
        Scalar objective = 0., max_violation = 0.;
        for (Index k = 0; k < info.dims.K; ++k)
        {
            const Index offset = info.offsets_tangent_u[k];
            const Index nu = info.dims.number_of_tangent_controls[k];
            for (Index j = 0; j < nu; ++j)
                objective += 0.1 * solution(offset + j) * solution(offset + j);
            for (Index i = 0; i < states; ++i)
            {
                const Scalar x = solution(offset + nu + i);
                objective += (x - 1.) * (x - 1.);
                if (k == 0) max_violation = std::max(max_violation, std::abs(x));
                if (k + 1 < info.dims.K)
                {
                    const Index next_offset = info.offsets_tangent_u[k + 1];
                    const Index next_nu = info.dims.number_of_tangent_controls[k + 1];
                    const Scalar u = solution(offset + i % 8);
                    const Scalar defect = 0.9 * x + 0.15 * u + 0.01 * u * u -
                                          solution(next_offset + next_nu + i);
                    max_violation = std::max(max_violation, std::abs(defect));
                }
            }
        }
        std::printf("states=%d mode=%s iterations=%d median_total_us=%.3f median_per_iteration_us=%.3f search_dir_us_per_iteration=%.3f hessian_eval_us_per_iteration=%.3f objective=%.9f max_eq_violation=%.3e\n",
                    int(states), mode, int(iterations), total_us[middle],
                    time_per_iteration_us[middle], 1e6 * timings.compute_search_dir.elapsed() / iterations,
                    1e6 * timings.eval_hessian.elapsed() / iterations, objective, max_violation);
    }
};
} // namespace

int main()
{
    for (Index states : {Index(16), Index(32), Index(64)})
    {
        Runner exact(states, "exact"), limited(states, "limited-memory"), full(states, "bfgs");
        for (int i = 0; i < 2; ++i) { exact.run(false); limited.run(false); full.run(false); }
        for (int i = 0; i < 11; ++i)
        {
            if (i & 1) { full.run(true); limited.run(true); exact.run(true); }
            else { exact.run(true); limited.run(true); full.run(true); }
        }
        exact.print(states, "exact");
        limited.print(states, "limited-memory");
        full.print(states, "bfgs");
    }
    return 0;
}
