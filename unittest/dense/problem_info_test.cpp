#include "fatrop/dense/problem_info.hpp"
#include "gtest/gtest.h"

namespace fatrop
{
namespace test
{

TEST(DenseProblemInfoTest, ConstructorAndOffsets)
{
    ProblemDims<DenseType> dims(/*nx=*/8, /*ng=*/3, /*ng_ineq=*/4);
    ProblemInfo<DenseType> info(dims);

    EXPECT_EQ(info.number_of_primal_variables, 8);
    EXPECT_EQ(info.number_of_tangent_variables, 8);

    EXPECT_EQ(info.number_of_slack_variables, 4);
    EXPECT_EQ(info.offset_slack, 8);

    EXPECT_EQ(info.number_of_g_eq_slack, 4);
    EXPECT_EQ(info.number_of_eq_constraints, 7);
    EXPECT_EQ(info.offset_g_eq_slack, 3);

    // pd-orig block offsets: [primal | slack | mult | zl | zu]
    EXPECT_EQ(info.pd_orig_offset_primal, 0);
    EXPECT_EQ(info.pd_orig_offset_slack, 8);
    EXPECT_EQ(info.pd_orig_offset_mult, 12);
    EXPECT_EQ(info.pd_orig_offset_zl, 19);
    EXPECT_EQ(info.pd_orig_offset_zu, 23);

    // restoration: extra n+p slacks per equality constraint
    EXPECT_EQ(info.number_of_slack_variables_resto, 4 + 2 * 7);

    EXPECT_EQ(info.constraint_allows_dual_damping.size(), 7u);
    for (bool v : info.constraint_allows_dual_damping)
        EXPECT_TRUE(v);
}

TEST(DenseProblemInfoTest, NoInequalities)
{
    ProblemDims<DenseType> dims(/*nx=*/5, /*ng=*/2, /*ng_ineq=*/0);
    ProblemInfo<DenseType> info(dims);
    EXPECT_EQ(info.number_of_slack_variables, 0);
    EXPECT_EQ(info.number_of_eq_constraints, 2);
    EXPECT_EQ(info.offset_g_eq_slack, 2); // empty block, offset only
}

} // namespace test
} // namespace fatrop
