#include "../../examples/limited_memory_problem.hpp"
#include <gtest/gtest.h>
#include "fatrop/dense/aug_system_solver.hpp"
#include "fatrop/dense/pd_solver_orig.hpp"
#include "fatrop/ip_algorithm/ip_eq_mult_initializer.hpp"
#include "fatrop/ip_algorithm/ip_initializer.hpp"
#include "fatrop/ip_algorithm/ip_nlp_resto.hpp"
#include "fatrop/ip_algorithm/limited_memory_hessian_evaluator.hpp"
#include "fatrop/ip_algorithm/hessian_approximation_evaluator.hpp"
#include "fatrop/ip_algorithm/limited_memory_nlp_resto.hpp"

namespace
{
using namespace limited_memory_example;

// Inspect extension state without adding approximation fields to native Hessians
// or algorithm data.
const LimitedMemoryBfgs<DenseType> &approximation(IpData<DenseType> &data)
{
    auto &evaluator = data.hessian_evaluator();
    if (auto *selector = dynamic_cast<HessianApproximationEvaluator<DenseType> *>(&evaluator))
        return selector->limited_memory();
    return dynamic_cast<LimitedMemoryHessianEvaluator<DenseType> &>(evaluator).approximation();
}
const std::shared_ptr<const LowRankHessian> &model(const Hessian<DenseType> &hessian)
{
    return dynamic_cast<const HessianWithLowRank<DenseType> &>(hessian).low_rank;
}
void check(const Result &r)
{
    EXPECT_TRUE(r.success);
    EXPECT_LT(r.max_error, 1e-6);
    EXPECT_LT(r.feasibility, 1e-8);
    EXPECT_LT(r.iterations, 100);
}

TEST(LimitedMemoryRegression, NonlinearEqualitiesMatchDenseAndOcp)
{
    auto dense = std::make_shared<DenseProblem>();
    auto ocp = std::make_shared<OcpProblem>();
    const auto d = solve_fatrop<DenseType>(std::make_shared<NlpDense>(dense), dense->model);
    const auto o = solve_fatrop<OcpType>(std::make_shared<NlpOcp>(ocp), ocp->dense.model);
    check(d); check(o);
    EXPECT_EQ(dense->hessian_calls, 0);
    EXPECT_EQ(ocp->dense.hessian_calls, 0);
    EXPECT_EQ(d.iterations, o.iterations);
    for (Index i = 0; i < Model::n; ++i) EXPECT_NEAR(d.x[i], o.x[i], 1e-9);
}

TEST(LimitedMemoryRegression, NonlinearGraphMatchesDense)
{
    auto dense = std::make_shared<DenseProblem>(false);
    auto graph = std::make_shared<limited_memory_example::GraphProblem>();
    const auto d = solve_fatrop<DenseType>(std::make_shared<NlpDense>(dense), dense->model);
    const auto g = solve_fatrop<fatrop::GraphProblem>(std::make_shared<NlpGraph>(graph), graph->model);
    check(d); check(g);
    EXPECT_EQ(dense->hessian_calls, 0);
    EXPECT_EQ(graph->hessian_calls, 0);
    EXPECT_EQ(d.iterations, g.iterations);
    for (Index i = 0; i < Model::n; ++i) EXPECT_NEAR(d.x[i], g.x[i], 1e-9);
}

TEST(LimitedMemoryRegression, FullBfgsDenseOcpAndSparseGraph)
{
    auto dense = std::make_shared<DenseProblem>();
    auto ocp = std::make_shared<OcpProblem>();
    auto graph = std::make_shared<limited_memory_example::GraphProblem>();
    const auto d = solve_fatrop<DenseType>(std::make_shared<NlpDense>(dense), dense->model, 6, "bfgs");
    const auto o = solve_fatrop<OcpType>(std::make_shared<NlpOcp>(ocp), ocp->dense.model, 6, "bfgs");
    const auto g = solve_fatrop<fatrop::GraphProblem>(std::make_shared<NlpGraph>(graph), graph->model, 6, "bfgs");
    check(d); check(o); check(g);
    EXPECT_EQ(dense->hessian_calls, 0);
    EXPECT_EQ(ocp->dense.hessian_calls, 0);
    EXPECT_EQ(graph->hessian_calls, 0);
}

TEST(LimitedMemoryRegression, OnePairHistoryConverges)
{
    auto dense = std::make_shared<DenseProblem>();
    const auto d = solve_fatrop<DenseType>(std::make_shared<NlpDense>(dense), dense->model, 1);
    check(d);
    EXPECT_EQ(dense->hessian_calls, 0);
}

class UnconstrainedDenseProblem : public DenseProblem
{
public:
    UnconstrainedDenseProblem() : DenseProblem(false) {}
    Index get_ng_ineq() const override { return 0; }
    Index eval_Ggt_ineq(const Scalar *, MAT *) override { return 0; }
    Index eval_gineq(const Scalar *, Scalar *) override { return 0; }
    Index get_bounds(Scalar *, Scalar *) const override { return 0; }
};

TEST(LimitedMemoryRegression, UnconstrainedProblemWithoutSlackOrMultiplierBlocks)
{
    auto problem = std::make_shared<UnconstrainedDenseProblem>();
    auto result = solve_fatrop<DenseType>(std::make_shared<NlpDense>(problem), problem->model);
    check(result);
    EXPECT_EQ(problem->hessian_calls, 0);
}

TEST(LimitedMemoryRegression, NonlinearMultistageDynamicsAndRepeatedOptimization)
{
    auto problem = std::make_shared<DynamicsProblem>();
    OptionRegistry options;
    IpAlgBuilder<OcpType> builder(std::make_shared<NlpOcp>(problem));
    auto algorithm = builder.with_options_registry(&options).build();
    options.set_option("hessian_approximation", std::string("limited-memory"));
    options.set_option("tolerance", Scalar(1e-8));
    options.set_option("constr_viol_tol", Scalar(1e-8));
    options.set_option("compl_inf_tol", Scalar(1e-8));
    options.set_option("tol_acceptable", Scalar(1e-8));
    options.set_option("max_iter", Index(150));
    options.set_option("print_level", Index(0));
    Index first_iterations = -1;
    for (Index run = 0; run < 2; ++run)
    {
        ASSERT_EQ(algorithm->optimize(), IpSolverReturnFlag::Success);
        auto data = builder.get_ipdata();
        const auto &x = data->current_iterate().primal_x();
        Index offset = 0;
        for (Index k = 0; k < DynamicsProblem::horizon; ++k)
        {
            if (problem->get_nu(k)) EXPECT_NEAR(x(offset++), 0, 1e-6);
            EXPECT_NEAR(x(offset++), 1, 1e-6);
        }
        EXPECT_EQ(problem->hessian_calls, 0);
        if (run == 0) first_iterations = data->iteration_number();
        else EXPECT_EQ(data->iteration_number(), first_iterations);
    }
}

TEST(LimitedMemoryRegression, FullBfgsMultistageOcp)
{
    auto problem = std::make_shared<DynamicsProblem>();
    OptionRegistry options;
    IpAlgBuilder<OcpType> builder(std::make_shared<NlpOcp>(problem));
    auto algorithm = builder.with_options_registry(&options).build();
    options.set_option("hessian_approximation", std::string("bfgs"));
    options.set_option("tolerance", Scalar(1e-8));
    options.set_option("max_iter", Index(150));
    options.set_option("print_level", Index(0));
    ASSERT_EQ(algorithm->optimize(), IpSolverReturnFlag::Success);
    const auto &x = builder.get_ipdata()->current_iterate().primal_x();
    Index offset = 0;
    for (Index k = 0; k < DynamicsProblem::horizon; ++k)
    {
        if (problem->get_nu(k)) EXPECT_NEAR(x(offset++), 0., 1e-5);
        EXPECT_NEAR(x(offset++), 1., 1e-5);
    }
    EXPECT_EQ(problem->hessian_calls, 0);
}

TEST(LimitedMemoryRegression, InitialDualsSolveLeastSquaresWithoutHessianOrHistory)
{
    auto problem = std::make_shared<DenseProblem>();
    OptionRegistry options;
    IpAlgBuilder<DenseType> builder(std::make_shared<NlpDense>(problem));
    auto algorithm = builder.with_options_registry(&options).build();
    options.set_option("hessian_approximation", std::string("limited-memory"));
    options.set_option("limited_memory_init_val", Scalar(3));
    auto data = builder.get_ipdata();
    auto dual_initializer = std::make_shared<IpEqMultInitializer<DenseType>>(data, builder.get_pd_solver());
    IpInitializer<DenseType> initializer(data, dual_initializer);
    initializer.initialize();
    EXPECT_EQ(problem->hessian_calls, 0);
    EXPECT_EQ(approximation(*data).history_size(), 0);
    EXPECT_EQ(approximation(*data).skipped_updates(), 0);

    // A * stationarity_residual = 0 are the normal equations for the
    // identity-on-(x,s) KKT multiplier least-squares problem.
    auto &it = data->current_iterate();
    const auto &rx = it.dual_infeasibility_x();
    const auto &rs = it.dual_infeasibility_s();
    const auto &x = it.primal_x();
    Scalar eq_residual = 0;
    for (Index i = 0; i < Model::n; ++i)
    {
        eq_residual += (i == 0 ? 2 * x(0) : 1) * rx(i);
        EXPECT_NEAR(Model::inequality_derivative(x(i)) * rx(i) - rs(i), 0, 1e-12);
        EXPECT_EQ(it.dual_bounds_l()(i), 1);
        EXPECT_EQ(it.dual_bounds_u()(i), 1);
    }
    EXPECT_NEAR(eq_residual, 0, 1e-12);
    const auto &initial_hessian = it.hessian();
    ASSERT_TRUE(model(initial_hessian));
    EXPECT_EQ(model(initial_hessian)->rank(), 0);
    EXPECT_EQ(model(initial_hessian)->sigma, 3);

    // A genuine displacement creates a pair; restarting discards it and
    // restores the configured initial scalar without evaluating exact H.
    VecRealAllocated next_x(Model::n);
    next_x = x;
    next_x(Model::n - 1) += 0.1;
    it.set_primal_x(next_x);
    it.hessian();
    EXPECT_EQ(approximation(*data).history_size(), 1);
    data->reset();
    initializer.initialize();
    EXPECT_EQ(approximation(*data).history_size(), 0);
    EXPECT_EQ(model(data->current_iterate().hessian())->sigma, 3);
    EXPECT_EQ(problem->hessian_calls, 0);
}

TEST(LimitedMemoryRegression, RestorationUsesBfgsWithoutDoubleCountingProximityCurvature)
{
    auto problem = std::make_shared<DenseProblem>();
    auto nlp = std::make_shared<NlpDense>(problem);
    auto evaluator = std::make_shared<LimitedMemoryHessianEvaluator<DenseType>>(ProblemInfo<DenseType>(nlp->problem_dims()));
    auto resto = std::make_shared<LimitedMemoryNlpResto<DenseType>>(nlp, evaluator);
    auto data = std::make_shared<IpData<DenseType>>(resto, evaluator);
    OptionRegistry options;
    data->register_options(options);
    resto->register_options(options);
    options.set_option("hessian_approximation", std::string("limited-memory"));
    data->reset(true);
    VecRealAllocated x(Model::n), s(Model::n);
    for (Index i = 0; i < Model::n; ++i) { x(i) = 0.5; s(i) = Model::inequality(x(i)); }
    resto->set_xs_reference(data->info(), x, s);
    resto->set_zeta(0.25);
    auto &it = data->current_iterate();
    it.set_primal_x(x);
    it.set_primal_s(VecRealScalar(it.primal_s().m(), 0.));
    it.set_dual_eq(VecRealScalar(it.dual_eq().m(), 0.));
    const auto &h = it.hessian();
    ASSERT_TRUE(model(h));
    EXPECT_EQ(model(h)->rank(), 0);
    EXPECT_EQ(model(h)->sigma, 1.);
    const auto &damping = it.primal_damping();
    for (Index i = 0; i < Model::n; ++i) EXPECT_EQ(damping(i), 0.);
    for (Index i = Model::n; i < 2 * Model::n; ++i) EXPECT_EQ(damping(i), 0.25);
    for (Index i = 2 * Model::n; i < damping.m(); ++i) EXPECT_EQ(damping(i), 0);
    EXPECT_EQ(approximation(*data).history_size(), 0);
    x(0) += .1;
    it.set_primal_x(x);
    const auto &updated = it.hessian();
    EXPECT_EQ(approximation(*data).history_size(), 1);
    EXPECT_NEAR(model(updated)->sigma, .25, 1e-14);
    VecRealAllocated step(Model::n), product(Model::n);
    step = 0.; step(0) = 1.; product = 0.;
    updated.apply_on_right(data->info(), step, 0., product, product);
    EXPECT_NEAR(product(0) + it.primal_damping()(0), .25, 1e-14);
    EXPECT_EQ(problem->hessian_calls, 0);
}

TEST(LimitedMemoryRegression, TrialEvaluationAndWatchdogRollbackPreserveHistory)
{
    auto problem = std::make_shared<DenseProblem>();
    OptionRegistry options;
    IpAlgBuilder<DenseType> builder(std::make_shared<NlpDense>(problem));
    auto algorithm = builder.with_options_registry(&options).build();
    options.set_option("hessian_approximation", std::string("limited-memory"));
    options.set_option("limited_memory_init_val", Scalar(3));
    auto data = builder.get_ipdata();
    auto dual_initializer = std::make_shared<IpEqMultInitializer<DenseType>>(data, builder.get_pd_solver());
    IpInitializer<DenseType> initializer(data, dual_initializer);
    initializer.initialize();
    auto initial = model(data->current_iterate().hessian());

    auto prepare_trial = [&](Scalar displacement) {
        auto &current = data->current_iterate();
        auto &trial = data->trial_iterate();
        VecRealAllocated x(Model::n);
        x = current.primal_x();
        x(Model::n - 1) += displacement;
        trial.set_primal_x(x);
        trial.set_primal_s(current.primal_s());
        trial.set_dual_eq(current.dual_eq());
        trial.set_dual_bounds_l(current.dual_bounds_l());
        trial.set_dual_bounds_u(current.dual_bounds_u());
    };
    prepare_trial(.1);
    EXPECT_EQ(model(data->trial_iterate().hessian()), initial);
    EXPECT_EQ(approximation(*data).history_size(), 0);
    EXPECT_EQ(approximation(*data).skipped_updates(), 0);
    data->accept_trial_iterate();
    auto saved = model(data->current_iterate().hessian());
    EXPECT_EQ(approximation(*data).history_size(), 1);

    // The watchdog stores the model matching its stored search direction.
    data->store_current_iterate();
    prepare_trial(.1);
    data->accept_trial_iterate();
    data->current_iterate().hessian();
    EXPECT_EQ(approximation(*data).history_size(), 2);
    data->restore_current_iterate();
    EXPECT_EQ(model(data->current_iterate().hessian()), saved);
    EXPECT_EQ(approximation(*data).history_size(), 0);

    // Its next accepted iterate starts a fresh sequence, even if the trial
    // Hessian was requested while retrying the restored direction.
    prepare_trial(.05);
    EXPECT_EQ(model(data->trial_iterate().hessian()), saved);
    data->accept_trial_iterate();
    auto restarted = model(data->current_iterate().hessian());
    EXPECT_DOUBLE_EQ(restarted->sigma, 3.);
    EXPECT_EQ(restarted->rank(), 0);
    EXPECT_EQ(problem->hessian_calls, 0);
}

class ExactDenseProblem : public DenseProblem
{
public:
    Index eval_Hh(const Scalar *scale, const Scalar *x, const Scalar *lambda, MAT *r) override
    {
        ++hessian_calls;
        blasfeo_gese_wrap(r->m, r->n, 0, r, 0, 0);
        for (Index i = 0; i < Model::n; ++i)
        {
            const Scalar e = x[i] - 1;
            Scalar curvature = *scale * (i + 1 + 3 * e * e) + 0.2 * lambda[1 + i];
            if (i == 0) curvature += 2 * (lambda[0] - *scale);
            blasfeo_matel_wrap(r, i, i) = curvature;
            blasfeo_matel_wrap(r, Model::n, i) = *scale * model.gradient(i, x);
        }
        return 0;
    }
};

TEST(LimitedMemoryRegression, SameSolverCanSwitchBackToExactHessian)
{
    auto problem = std::make_shared<ExactDenseProblem>();
    OptionRegistry options;
    IpAlgBuilder<DenseType> builder(std::make_shared<NlpDense>(problem));
    auto algorithm = builder.with_options_registry(&options).build();
    options.set_option("tolerance", Scalar(1e-8));
    options.set_option("print_level", Index(0));
    options.set_option("hessian_approximation", std::string("limited-memory"));
    ASSERT_EQ(algorithm->optimize(), IpSolverReturnFlag::Success);
    EXPECT_EQ(problem->hessian_calls, 0);
    EXPECT_GT(approximation(*builder.get_ipdata()).history_size(), 0);
    options.set_option("hessian_approximation", std::string("bfgs"));
    ASSERT_EQ(algorithm->optimize(), IpSolverReturnFlag::Success);
    EXPECT_EQ(problem->hessian_calls, 0);
    EXPECT_EQ(approximation(*builder.get_ipdata()).history_size(), 0);
    options.set_option("hessian_approximation", std::string("exact"));
    ASSERT_EQ(algorithm->optimize(), IpSolverReturnFlag::Success);
    EXPECT_GT(problem->hessian_calls, 0);
    const auto &x = builder.get_ipdata()->current_iterate().primal_x();
    for (Index i = 0; i < Model::n; ++i) EXPECT_NEAR(x(i), 1, 1e-6);
    EXPECT_EQ(approximation(*builder.get_ipdata()).history_size(), 0);
}
TEST(LimitedMemoryRegression, ConstructionCapacitySupportsLargerHistoryAndNarrowRhsPanels)
{
    auto problem = std::make_shared<DenseProblem>();
    IpWorkspace workspace;
    workspace.hessian_history_capacity = 9;
    workspace.rhs_capacity = 3;
    OptionRegistry options;
    IpAlgBuilder<DenseType> builder(std::make_shared<NlpDense>(problem), workspace);
    auto algorithm = builder.with_options_registry(&options).build();
    options.set_option("hessian_approximation", std::string("limited-memory"));
    options.set_option("limited_memory_max_history", Index(9));
    options.set_option("print_level", Index(0));
    ASSERT_EQ(algorithm->optimize(), IpSolverReturnFlag::Success);
    EXPECT_EQ(approximation(*builder.get_ipdata()).history_size(), 9);
    const auto &x = builder.get_ipdata()->current_iterate().primal_x();
    for (Index i = 0; i < Model::n; ++i) EXPECT_NEAR(x(i), 1., 1e-6);
    EXPECT_EQ(problem->hessian_calls, 0);
    EXPECT_THROW(options.set_option("limited_memory_max_history", Index(10)), std::runtime_error);
}
} // namespace
