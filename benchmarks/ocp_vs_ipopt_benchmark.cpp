#define main fatrop_hessian_mode_benchmark_main
#include "ocp_hessian_mode_benchmark.cpp"
#undef main

#include <IpIpoptApplication.hpp>
#include <IpSolveStatistics.hpp>
#include <IpTNLP.hpp>
#include <numeric>

namespace
{
class IpoptOcp final : public Ipopt::TNLP
{
public:
    IpoptOcp(int states, int stages = 60) : states_(states), stages_(stages) {}
    int iterations = 0;
    bool success = false;
    double objective = 0.;
    double max_eq_violation = 0.;

    int variables() const { return stages_ * states_ + (stages_ - 1) * 8; }
    int equalities() const { return stages_ * states_; }
    int controls() const { return (stages_ - 1) * 8; }
    int stage_offset(int k) const { return k * (states_ + 8); }
    int state_offset(int k) const { return stage_offset(k) + (k + 1 < stages_ ? 8 : 0); }
    int initial_row(int i) const { return i; }
    int dynamic_row(int k, int i) const { return states_ + k * states_ + i; }
    int control_row(int k, int j) const { return equalities() + 8 * k + j; }

    bool get_nlp_info(Ipopt::Index &n, Ipopt::Index &m, Ipopt::Index &nnz_jac,
                      Ipopt::Index &nnz_h, IndexStyleEnum &style) override
    {
        n = variables(); m = equalities() + controls();
        nnz_jac = states_ + (stages_ - 1) * states_ * 3 + controls();
        nnz_h = 0; style = C_STYLE; return true;
    }
    bool get_bounds_info(Ipopt::Index n, Ipopt::Number *xl, Ipopt::Number *xu,
                         Ipopt::Index m, Ipopt::Number *gl, Ipopt::Number *gu) override
    {
        for (int i = 0; i < n; ++i) { xl[i] = -1e20; xu[i] = 1e20; }
        for (int i = 0; i < equalities(); ++i) gl[i] = gu[i] = 0.;
        for (int i = equalities(); i < m; ++i) { gl[i] = -2.; gu[i] = 2.; }
        return true;
    }
    bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number *x,
                            bool init_z, Ipopt::Number *, Ipopt::Number *,
                            Ipopt::Index, bool init_lambda, Ipopt::Number *) override
    {
        if (!init_x || init_z || init_lambda) return false;
        std::fill(x, x + n, 0.); return true;
    }
    bool eval_f(Ipopt::Index, const Ipopt::Number *x, bool, Ipopt::Number &f) override
    {
        f = 0.;
        for (int k = 0; k < stages_; ++k)
        {
            if (k + 1 < stages_)
                for (int j = 0; j < 8; ++j) f += 0.1 * x[stage_offset(k) + j] * x[stage_offset(k) + j];
            for (int i = 0; i < states_; ++i)
            {
                const double e = x[state_offset(k) + i] - 1.; f += e * e;
            }
        }
        return true;
    }
    bool eval_grad_f(Ipopt::Index, const Ipopt::Number *x, bool, Ipopt::Number *grad) override
    {
        for (int k = 0; k < stages_; ++k)
        {
            if (k + 1 < stages_)
                for (int j = 0; j < 8; ++j) grad[stage_offset(k) + j] = 0.2 * x[stage_offset(k) + j];
            for (int i = 0; i < states_; ++i)
                grad[state_offset(k) + i] = 2. * (x[state_offset(k) + i] - 1.);
        }
        return true;
    }
    bool eval_g(Ipopt::Index, const Ipopt::Number *x, bool, Ipopt::Index,
                Ipopt::Number *g) override
    {
        for (int i = 0; i < states_; ++i) g[initial_row(i)] = x[state_offset(0) + i];
        for (int k = 0; k + 1 < stages_; ++k)
        {
            for (int i = 0; i < states_; ++i)
            {
                const double u = x[stage_offset(k) + i % 8];
                g[dynamic_row(k, i)] = 0.9 * x[state_offset(k) + i] + 0.15 * u
                                    + 0.01 * u * u - x[state_offset(k + 1) + i];
            }
            for (int j = 0; j < 8; ++j) g[control_row(k, j)] = x[stage_offset(k) + j];
        }
        return true;
    }
    bool eval_jac_g(Ipopt::Index, const Ipopt::Number *x, bool, Ipopt::Index,
                    Ipopt::Index, Ipopt::Index *rows, Ipopt::Index *cols,
                    Ipopt::Number *values) override
    {
        int p = 0;
        auto entry = [&](int row, int col, double value) {
            if (values) values[p] = value;
            else { rows[p] = row; cols[p] = col; }
            ++p;
        };
        for (int i = 0; i < states_; ++i) entry(initial_row(i), state_offset(0) + i, 1.);
        for (int k = 0; k + 1 < stages_; ++k)
        {
            for (int i = 0; i < states_; ++i)
            {
                const int j = i % 8;
                entry(dynamic_row(k, i), stage_offset(k) + j,
                      values ? 0.15 + 0.02 * x[stage_offset(k) + j] : 0.);
                entry(dynamic_row(k, i), state_offset(k) + i, 0.9);
                entry(dynamic_row(k, i), state_offset(k + 1) + i, -1.);
            }
            for (int j = 0; j < 8; ++j) entry(control_row(k, j), stage_offset(k) + j, 1.);
        }
        return true;
    }
    bool eval_h(Ipopt::Index, const Ipopt::Number *, bool, Ipopt::Number,
                Ipopt::Index, const Ipopt::Number *, bool, Ipopt::Index,
                Ipopt::Index *, Ipopt::Index *, Ipopt::Number *) override { return false; }
    void finalize_solution(Ipopt::SolverReturn status, Ipopt::Index,
                           const Ipopt::Number *, const Ipopt::Number *,
                           const Ipopt::Number *, Ipopt::Index m,
                           const Ipopt::Number *g, const Ipopt::Number *,
                           Ipopt::Number f, const Ipopt::IpoptData *,
                           Ipopt::IpoptCalculatedQuantities *) override
    {
        success = status == Ipopt::SUCCESS;
        objective = f;
        max_eq_violation = 0.;
        for (int i = 0; i < equalities(); ++i)
            max_eq_violation = std::max(max_eq_violation, std::abs(g[i]));
    }
private:
    int states_, stages_;
};

void compare(int states)
{
    Runner fatrop(states, "limited-memory");
    Runner exact(states, "exact");
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
    auto opts = app->Options();
    opts->SetStringValue("hessian_approximation", "limited-memory");
    opts->SetStringValue("limited_memory_initialization", "scalar1");
    opts->SetIntegerValue("limited_memory_max_history", 6);
    opts->SetStringValue("nlp_scaling_method", "none");
    opts->SetStringValue("mu_strategy", "monotone");
    opts->SetNumericValue("tol", 1e-8);
    opts->SetNumericValue("constr_viol_tol", 1e-8);
    opts->SetNumericValue("dual_inf_tol", 1e-8);
    opts->SetNumericValue("compl_inf_tol", 1e-8);
    opts->SetIntegerValue("acceptable_iter", 0);
    opts->SetIntegerValue("max_iter", 300);
    opts->SetIntegerValue("print_level", 0);
    opts->SetStringValue("sb", "yes");
    if (app->Initialize() != Ipopt::Solve_Succeeded) std::abort();
    Ipopt::SmartPtr<IpoptOcp> problem = new IpoptOcp(states);
    std::vector<double> ipopt_us;
    for (int trial = 0; trial < 7; ++trial)
    {
        fatrop.run(trial >= 2);
        exact.run(trial >= 2);
        const auto start = std::chrono::steady_clock::now();
        app->OptimizeTNLP(problem);
        const double us = std::chrono::duration<double, std::micro>(
            std::chrono::steady_clock::now() - start).count();
        if (!problem->success) std::abort();
        if (trial >= 2) ipopt_us.push_back(us);
    }
    std::sort(ipopt_us.begin(), ipopt_us.end());
    std::sort(fatrop.total_us.begin(), fatrop.total_us.end());
    std::sort(exact.total_us.begin(), exact.total_us.end());
    const double ft = fatrop.total_us[fatrop.total_us.size() / 2];
    const double et = exact.total_us[exact.total_us.size() / 2];
    const double it = ipopt_us[ipopt_us.size() / 2];
    std::printf("states=%d stages=60 variables=%d constraints=%d "
                "fatrop_bfgs_iter=%d fatrop_exact_iter=%d ipopt_bfgs_iter=%d "
                "fatrop_bfgs_median_ms=%.3f fatrop_exact_median_ms=%.3f "
                "ipopt_bfgs_median_ms=%.3f "
                "ratio_bfgs_over_exact=%.2f ratio_ipopt_over_exact=%.2f "
                "fatrop_bfgs_obj=%.8g fatrop_exact_obj=%.8g ipopt_obj=%.8g "
                "ipopt_eq_violation=%.2e\n",
                states, problem->variables(), problem->equalities() + problem->controls(),
                int(fatrop.iterations), int(exact.iterations),
                int(app->Statistics()->IterationCount()),
                ft / 1000., et / 1000., it / 1000., ft / et, it / et,
                fatrop.builder.get_ipdata()->current_iterate().obj_value(),
                exact.builder.get_ipdata()->current_iterate().obj_value(),
                problem->objective, problem->max_eq_violation);
}
}

int main()
{
    for (int states : {16, 64}) compare(states);
}
