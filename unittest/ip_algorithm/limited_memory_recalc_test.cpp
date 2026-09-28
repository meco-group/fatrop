#include "../../examples/limited_memory_problem.hpp"
#include <gtest/gtest.h>

namespace
{
using namespace limited_memory_example;

Index solve_chain_with_recalculation(bool set_option, bool enabled)
{
    auto problem = std::make_shared<DenseProblem>(false);
    OptionRegistry options;
    IpAlgBuilder<DenseType> builder(std::make_shared<NlpDense>(problem));
    auto algorithm = builder.with_options_registry(&options).build();
    options.set_option("hessian_approximation", std::string("limited-memory"));
    options.set_option("limited_memory_max_history", Index(6));
    options.set_option("tolerance", Scalar(1e-8));
    options.set_option("constr_viol_tol", Scalar(1e-8));
    options.set_option("compl_inf_tol", Scalar(1e-8));
    options.set_option("tol_acceptable", Scalar(1e-8));
    options.set_option("max_iter", Index(300));
    options.set_option("print_level", Index(0));
    if (set_option) options.set_option("recalc_y", enabled);
    EXPECT_EQ(algorithm->optimize(), IpSolverReturnFlag::Success);
    EXPECT_EQ(problem->hessian_calls, 0);
    return builder.get_ipdata()->iteration_number();
}

TEST(LimitedMemoryRecalculation, DefaultMatchesExplicitYesAndChangesTrajectory)
{
    const Index by_default = solve_chain_with_recalculation(false, false);
    const Index explicit_yes = solve_chain_with_recalculation(true, true);
    const Index explicit_no = solve_chain_with_recalculation(true, false);
    EXPECT_EQ(by_default, 14);
    EXPECT_EQ(explicit_yes, by_default);
    EXPECT_EQ(explicit_no, 18);
}
}
