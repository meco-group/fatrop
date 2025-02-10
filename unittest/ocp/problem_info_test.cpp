#include "fatrop/ocp/problem_info.hpp"
#include "gtest/gtest.h"

namespace fatrop
{
    namespace test
    {

class ProblemInfoTest : public ::testing::Test
{
protected:
    ProblemDims<OcpType> createSampleDims()
    {
        return ProblemDims<OcpType>(
            3, // 3 stages
            std::vector<Index>{2, 3, 4}, // number_of_controls
            std::vector<Index>{4, 5, 6}, // number_of_states
            std::vector<Index>{1, 2, 3}, // number_of_eq_constraints
            std::vector<Index>{6, 7, 8}  // number_of_ineq_constraints
        );
    }
};

TEST_F(ProblemInfoTest, ConstructorAndOffsets)
{
    ProblemDims<OcpType> dims = createSampleDims();
            ProblemInfo<OcpType> problem_info(dims);

            // Test dimensions
            EXPECT_EQ(problem_info.dims.K, 3);

            // Test primal offsets
            EXPECT_EQ(problem_info.offset_primal, 0);
            EXPECT_EQ(problem_info.offsets_primal_u, std::vector<Index>({0, 6, 14}));
            EXPECT_EQ(problem_info.offsets_primal_x, std::vector<Index>({2, 9, 18}));

            // // Test slack offsets
            EXPECT_EQ(problem_info.offset_slack, 24);
            EXPECT_EQ(problem_info.offsets_slack, std::vector<Index>({0, 6, 13}));

            // // Test equality constraint numbers
            EXPECT_EQ(problem_info.number_of_g_eq_dyn, 11);
            EXPECT_EQ(problem_info.number_of_g_eq_path, 6);
            EXPECT_EQ(problem_info.number_of_g_eq_slack, 21);

            // // Test equality constraint offsets
            EXPECT_EQ(problem_info.offset_g_eq_dyn, 6);
            EXPECT_EQ(problem_info.offset_g_eq_path, 0);
            EXPECT_EQ(problem_info.offset_g_eq_slack, 17);

            EXPECT_EQ(problem_info.offsets_g_eq_path, std::vector<Index>({0, 1, 3}));
            EXPECT_EQ(problem_info.offsets_g_eq_dyn, std::vector<Index>({6, 11}));
            EXPECT_EQ(problem_info.offsets_g_eq_slack, std::vector<Index>({17, 23, 30}));
        }

    } // namespace test
} // namespace fatrop
