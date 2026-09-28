#include "../allocation_tracker.hpp"
#include "../../examples/limited_memory_problem.hpp"
#include <gtest/gtest.h>
#include <cstdlib>
#include <memory>
#include "fatrop/ip_algorithm/limited_memory_nlp_resto.hpp"

namespace
{
using namespace limited_memory_example;
namespace tracker = fatrop::test::allocation_tracker;

TEST(AllocationTracker, ObservesCppAndBlasfeoAllocations)
{
    // Sanity-check both C++ new and BLASFEO's C allocator. The latter allocation
    // originates in libfatrop / libblasfeo, proving shared-library visibility.
    tracker::start();
    auto *object = new VecRealAllocated(37);
    const auto counts = tracker::stop();
    delete object;
    EXPECT_GE(counts.calls, 2u);
    EXPECT_GE(counts.bytes, 37u * sizeof(Scalar));
}

template <typename ProblemType>
void check_optimization_allocations(std::shared_ptr<Nlp<ProblemType>> nlp,
                                    const IpWorkspace &workspace = {}, Index history = 6)
{
    OptionRegistry options;
    IpAlgBuilder<ProblemType> builder(nlp, workspace);
    auto algorithm = builder.with_options_registry(&options).build();
    options.set_option("hessian_approximation", std::string("limited-memory"));
    options.set_option("limited_memory_max_history", history);
    options.set_option("tolerance", Scalar(1e-8));
    options.set_option("constr_viol_tol", Scalar(1e-8));
    options.set_option("print_level", Index(0));
    options.set_option("suppress_banner", true);
    const bool trace = std::getenv("FATROP_ALLOC_TRACE") != nullptr;
    for (int run = 0; run < 2; ++run)
    {
        // This includes reset(), initial dual computation, the first KKT
        // factorization, every BFGS history growth/eviction and solve/refinement.
        tracker::start(trace);
        const auto status = algorithm->optimize();
        const auto counts = tracker::stop();
        EXPECT_EQ(status, IpSolverReturnFlag::Success);
        EXPECT_EQ(counts.calls, 0u) << "run " << run << ", bytes " << counts.bytes;
    }
}

TEST(LimitedMemoryAllocations, DenseFirstAndRepeatedOptimization)
{
    auto problem = std::make_shared<DenseProblem>();
    check_optimization_allocations<DenseType>(std::make_shared<NlpDense>(problem));
    EXPECT_EQ(problem->hessian_calls, 0);
}
TEST(LimitedMemoryAllocations, OcpFirstAndRepeatedOptimization)
{
    auto problem = std::make_shared<DynamicsProblem>();
    check_optimization_allocations<OcpType>(std::make_shared<NlpOcp>(problem));
    EXPECT_EQ(problem->hessian_calls, 0);
}
TEST(LimitedMemoryAllocations, GraphFirstAndRepeatedOptimization)
{
    auto problem = std::make_shared<limited_memory_example::GraphProblem>();
    check_optimization_allocations<fatrop::GraphProblem>(std::make_shared<NlpGraph>(problem));
    EXPECT_EQ(problem->hessian_calls, 0);
}
TEST(LimitedMemoryAllocations, ConfiguredLargerHistoryAndSmallRhsPanel)
{
    auto problem = std::make_shared<DenseProblem>();
    IpWorkspace workspace;
    workspace.hessian_history_capacity = 10;
    workspace.rhs_capacity = 3;
    workspace.max_iterations_capacity = 200;
    check_optimization_allocations<DenseType>(std::make_shared<NlpDense>(problem), workspace, 9);
    EXPECT_EQ(problem->hessian_calls, 0);
}

TEST(LimitedMemoryAllocations, SkippedUpdatesAndResetsReuseHistoryStorage)
{
    ProblemDims<DenseType> dims(1, 0, 0);
    ProblemInfo<DenseType> info(dims);
    Jacobian<DenseType> jac(dims);
    LimitedMemoryBfgs<DenseType> history(info);
    VecRealAllocated x(1), g(1), lambda(0);
    std::shared_ptr<const LowRankHessian> last;
    Index skipped = 0;
    tracker::start();
    for (int run = 0; run < 3; ++run)
    {
        history.reset();
        x = 0.; g = 0.; last = history.update(info, x, g, jac, lambda);
        x = 1.; g = 3.; last = history.update(info, x, g, jac, lambda);
        x = 2.; g = 2.; last = history.update(info, x, g, jac, lambda);
        x = 3.; g = 1.; last = history.update(info, x, g, jac, lambda);
        skipped = history.skipped_updates();
        x = 4.; g = 8.; last = history.update(info, x, g, jac, lambda);
    }
    const auto counts = tracker::stop();
    EXPECT_EQ(counts.calls, 0u);
    EXPECT_EQ(skipped, 2);
    EXPECT_EQ(history.history_size(), 0);
    EXPECT_EQ(last->rank(), 0);
}

TEST(LimitedMemoryAllocations, RestorationHistoryAndProximityReuseStorage)
{
    auto problem = std::make_shared<DenseProblem>();
    auto nlp = std::make_shared<NlpDense>(problem);
    auto evaluator = std::make_shared<LimitedMemoryHessianEvaluator<DenseType>>(
        ProblemInfo<DenseType>(nlp->problem_dims()));
    auto resto = std::make_shared<LimitedMemoryNlpResto<DenseType>>(nlp, evaluator);
    auto data = std::make_shared<IpData<DenseType>>(resto, evaluator);
    OptionRegistry options;
    data->register_options(options);
    options.set_option("hessian_approximation", std::string("limited-memory"));
    VecRealAllocated x(Model::n), s(Model::n);
    x = .5; s = Model::inequality(.5);
    resto->set_xs_reference(data->info(), x, s);
    resto->set_zeta(.25);
    tracker::start();
    for (int run = 0; run < 2; ++run)
    {
        data->reset(true);
        auto &it = data->current_iterate();
        it.set_primal_x(x);
        it.set_primal_s(VecRealScalar(it.primal_s().m(), 0.));
        it.set_dual_eq(VecRealScalar(it.dual_eq().m(), 0.));
        it.hessian();
        it.primal_damping();
        x(0) += .1;
        it.set_primal_x(x);
        it.hessian();
        it.primal_damping();
    }
    const auto counts = tracker::stop();
    EXPECT_EQ(counts.calls, 0u);
    EXPECT_EQ(evaluator->approximation().history_size(), 1);
    EXPECT_EQ(problem->hessian_calls, 0);
}

TEST(LimitedMemoryAllocations, IterationCapacityAndRuntimeOptions)
{
    auto problem = std::make_shared<DenseProblem>();
    IpWorkspace workspace;
    workspace.max_iterations_capacity = 60;
    OptionRegistry options;
    IpAlgBuilder<DenseType> builder(std::make_shared<NlpDense>(problem), workspace);
    auto algorithm = builder.with_options_registry(&options).build();
    options.set_option("hessian_approximation", std::string("limited-memory"));
    options.set_option("print_level", Index(0));
    options.set_option("suppress_banner", true);
    const std::string history_option = "limited_memory_max_history";
    tracker::start();
    options.set_option("max_iter", Index(60));
    options.set_option(history_option, Index(4));
    const auto status = algorithm->optimize();
    const auto counts = tracker::stop();
    EXPECT_EQ(status, IpSolverReturnFlag::Success);
    EXPECT_EQ(counts.calls, 0u);
    EXPECT_THROW(options.set_option("max_iter", Index(61)), std::runtime_error);
}
} // namespace
