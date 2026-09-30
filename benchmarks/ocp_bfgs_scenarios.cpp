#include "ocp_derivative_checker.hpp"
#include "problems/timed_problem.hpp"
#include "problems/problem_catalog.hpp"
#include "compiled_derivatives.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

namespace
{
using namespace fatrop;
using namespace fatrop_benchmark;
template<class P> using EvaluationProblem = CompiledProblem<P>;

struct Result
{
    Index iterations = 0;
    IpSolverReturnFlag status = IpSolverReturnFlag::Success;
    Scalar objective = 0.0, max_violation = 0.0;
    double total_ms = 0.0, evaluation_ms = 0.0, fatrop_ms = 0.0;
    std::int64_t evaluation_calls = 0;
    double hessian_ms = 0.0, dynamics_jacobian_ms = 0.0;
    std::int64_t hessian_calls = 0, dynamics_jacobian_calls = 0;
    std::array<double, 9> callback_ms{};
    std::array<std::int64_t, 9> callback_calls{};
};

template <class Problem> Result run(const char *mode, Index repeats)
{
    auto problem = std::make_shared<TimedProblem<EvaluationProblem<Problem>>>();
    OptionRegistry options;
    IpWorkspace workspace;
    if constexpr (std::is_same_v<Problem, MinimumTimeBicycleObstacle>)
        workspace.max_iterations_capacity = 3000;
    IpAlgBuilder<OcpType> builder(std::make_shared<NlpOcp>(problem), workspace);
    auto algorithm = builder.with_options_registry(&options).build();
    options.set_option("hessian_approximation", std::string(mode));
    options.set_option("tolerance", Scalar(1e-8));
    options.set_option("constr_viol_tol", Scalar(1e-8));
    options.set_option("compl_inf_tol", Scalar(1e-8));
    options.set_option("tol_acceptable", Scalar(1e-8));
    options.set_option("max_iter", Index(workspace.max_iterations_capacity));
    options.set_option("print_level", Index(0));
    options.set_option("suppress_banner", true);

    struct Sample
    {
        double total, evaluation, hessian, dynamics_jacobian;
        std::int64_t calls, hessian_calls, dynamics_jacobian_calls;
        std::array<double, 9> callback_ms{};
        std::array<std::int64_t, 9> callback_calls{};
    };
    std::vector<Sample> samples;
    Result result;
    Index reference_iterations = -1;
    for (Index run_index = 0; run_index <= repeats; ++run_index)
    {
        problem->reset_timing();
        const auto start = std::chrono::steady_clock::now();
        const auto status = algorithm->optimize();
        const double total_ms = std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - start).count();
        const double evaluation_ms = problem->evaluation_ms();
        const auto evaluation_calls = problem->evaluation_calls();
        result.status = status;
        result.iterations = builder.get_ipdata()->iteration_number();
        if (status != IpSolverReturnFlag::Success)
        {
            result.total_ms = total_ms;
            result.evaluation_ms = evaluation_ms;
            result.fatrop_ms = std::max(0.0, total_ms - evaluation_ms);
            result.evaluation_calls = evaluation_calls;
            result.hessian_ms = problem->callback_ms(TimedProblem<EvaluationProblem<Problem>>::Hessian);
            result.dynamics_jacobian_ms = problem->callback_ms(TimedProblem<EvaluationProblem<Problem>>::DynamicsJacobian);
            result.hessian_calls = problem->callback_calls(TimedProblem<EvaluationProblem<Problem>>::Hessian);
            result.dynamics_jacobian_calls = problem->callback_calls(TimedProblem<EvaluationProblem<Problem>>::DynamicsJacobian);
            for (int i = 0; i < 9; ++i)
            {
                const auto callback = static_cast<typename TimedProblem<EvaluationProblem<Problem>>::Callback>(i);
                result.callback_ms[i] = problem->callback_ms(callback);
                result.callback_calls[i] = problem->callback_calls(callback);
            }
            break;
        }
        if (run_index && reference_iterations != result.iterations)
        {
            std::fprintf(stderr, "iteration count changed between repeated solves\n");
            std::exit(1);
        }
        if (!run_index) reference_iterations = result.iterations;
        if (run_index)
        {
            Sample sample{total_ms, evaluation_ms,
                problem->callback_ms(TimedProblem<EvaluationProblem<Problem>>::Hessian),
                problem->callback_ms(TimedProblem<EvaluationProblem<Problem>>::DynamicsJacobian),
                evaluation_calls,
                problem->callback_calls(TimedProblem<EvaluationProblem<Problem>>::Hessian),
                problem->callback_calls(TimedProblem<EvaluationProblem<Problem>>::DynamicsJacobian)};
            for (int i = 0; i < 9; ++i)
            {
                const auto callback = static_cast<typename TimedProblem<EvaluationProblem<Problem>>::Callback>(i);
                sample.callback_ms[i] = problem->callback_ms(callback);
                sample.callback_calls[i] = problem->callback_calls(callback);
            }
            samples.push_back(sample);
        }
    }
    if (result.status == IpSolverReturnFlag::Success)
    {
        std::sort(samples.begin(), samples.end(),
            [](const Sample &a, const Sample &b) { return a.total < b.total; });
        const auto &middle = samples[samples.size() / 2];
        result.total_ms = middle.total;
        result.evaluation_ms = middle.evaluation;
        result.fatrop_ms = std::max(0.0, middle.total - middle.evaluation);
        result.evaluation_calls = middle.calls;
        result.hessian_ms = middle.hessian;
        result.dynamics_jacobian_ms = middle.dynamics_jacobian;
        result.hessian_calls = middle.hessian_calls;
        result.dynamics_jacobian_calls = middle.dynamics_jacobian_calls;
        result.callback_ms = middle.callback_ms;
        result.callback_calls = middle.callback_calls;
    }

    const auto data = builder.get_ipdata();
    const auto &info = data->info();
    const auto &solution = data->current_iterate().primal_x();
    for (Index k = 0; k < Problem::stages; ++k)
    {
        const Index nu = problem->get_nu(k);
        const Scalar *u = solution.data() + info.offsets_tangent_u[k];
        const Scalar *x = u + nu;
        result.objective += problem->stage_cost(u, x, k);
        if (k == 0)
            for (Index i = 0; i < Problem::states; ++i)
                result.max_violation = std::max(result.max_violation,
                    std::abs(x[i] - problem->initial(i)));
        if (k + 1 == Problem::stages)
            for (Index i = 0; i < Problem::terminal_equalities; ++i)
                result.max_violation = std::max(result.max_violation,
                    std::abs(x[i] - problem->target(i)));
        for (Index j = 0; j < nu; ++j)
            result.max_violation = std::max(result.max_violation,
                std::max(problem->lower_control(j) - u[j], u[j] - problem->upper_control(j)));
        for (Index q = 0; q < problem->extra_inequalities(k); ++q)
        {
            const Scalar value = problem->inequality_value(q, u, x, k);
            result.max_violation = std::max(result.max_violation,
                std::max(problem->lower_inequality(q, k) - value,
                         value - problem->upper_inequality(q, k)));
        }
        if (k + 1 < Problem::stages)
        {
            const Scalar *next = solution.data() + info.offsets_tangent_u[k + 1]
                               + problem->get_nu(k + 1);
            for (Index i = 0; i < Problem::states; ++i)
                result.max_violation = std::max(result.max_violation,
                    std::abs(problem->dynamics_value(i, u, x, k) - next[i]));
        }
    }
    return result;
}

template <class Problem> bool benchmark(Index repeats, const char *filter)
{
    if (filter && std::string(filter) != Problem::name) return true;
    EvaluationProblem<Problem> derivative_problem;
    OcpDerivativeChecker checker(derivative_problem);
    const auto report = checker.check_all_initial_points();
    std::printf("scenario=%s derivative_checks=%d derivative_failures=%d max_scaled_error=%.3e\n",
        Problem::name, report.checks, report.failures, report.max_scaled_error);
    if (report.failures)
    {
        std::fprintf(stderr, "worst derivative: %s\n", report.worst.c_str());
        return false;
    }
    for (const char *mode : {"exact", "bfgs", "limited-memory"})
    {
        const Result result = run<Problem>(mode, repeats);
        std::printf("scenario=%s mode=%s status=%d iterations=%d total_ms=%.3f evaluation_ms=%.3f fatrop_ms=%.3f evaluation_calls=%lld hessian_ms=%.3f hessian_calls=%lld dynamics_jacobian_ms=%.3f dynamics_jacobian_calls=%lld cost_gradient_ms=%.3f cost_gradient_calls=%lld inequality_jacobian_ms=%.3f inequality_jacobian_calls=%lld dynamics_value_ms=%.3f cost_value_ms=%.3f inequality_value_ms=%.3f objective=%.9f max_violation=%.3e\n",
            Problem::name, mode, int(result.status), int(result.iterations),
            result.total_ms, result.evaluation_ms, result.fatrop_ms,
            static_cast<long long>(result.evaluation_calls),
            result.hessian_ms, static_cast<long long>(result.hessian_calls),
            result.dynamics_jacobian_ms, static_cast<long long>(result.dynamics_jacobian_calls),
            result.callback_ms[TimedProblem<EvaluationProblem<Problem>>::CostGradient],
            static_cast<long long>(result.callback_calls[TimedProblem<EvaluationProblem<Problem>>::CostGradient]),
            result.callback_ms[TimedProblem<EvaluationProblem<Problem>>::InequalityJacobian],
            static_cast<long long>(result.callback_calls[TimedProblem<EvaluationProblem<Problem>>::InequalityJacobian]),
            result.callback_ms[TimedProblem<EvaluationProblem<Problem>>::DynamicsValue],
            result.callback_ms[TimedProblem<EvaluationProblem<Problem>>::CostValue],
            result.callback_ms[TimedProblem<EvaluationProblem<Problem>>::InequalityValue],
            result.objective, result.max_violation);
    }
    return true;
}
}

int main(int argc, char **argv)
{
    std::printf("ad_backend=compiled\n");
    const fatrop::Index repeats = argc > 1
        ? std::max(fatrop::Index(1), fatrop::Index(std::atoi(argv[1]))) : 7;
    const char *filter = argc > 2 ? argv[2] : nullptr;
#define BENCHMARK(Problem) if (!benchmark<fatrop_benchmark::Problem>(repeats, filter)) return 1;
    FATROP_BENCHMARK_PROBLEMS(BENCHMARK)
#undef BENCHMARK
}
