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
            EXPECT_EQ(problem_info.offsets_g_eq_dyn, std::vector<Index>({6, 11, 0}));
            EXPECT_EQ(problem_info.offsets_g_eq_slack, std::vector<Index>({17, 23, 30}));
        }

TEST_F(ProblemInfoTest, ManifoldSlackOffsetsAreTangentSized)
{
    ProblemDims<OcpType> dims(
        2,
        std::vector<Index>{2, 0},
        std::vector<Index>{4, 5},
        std::vector<Index>{0, 0},
        std::vector<Index>{3, 0},
        std::vector<Index>{2, 0},
        std::vector<Index>{3, 4});
    ProblemInfo<OcpType> info(dims);

    const Index n_primal = 4 + 2 + 5;
    const Index n_tangent = 3 + 2 + 4;
    const Index n_slack = 3;
    EXPECT_EQ(info.number_of_primal_variables, n_primal);
    EXPECT_EQ(info.number_of_tangent_variables, n_tangent);
    EXPECT_EQ(info.number_of_slack_variables, n_slack);
    EXPECT_EQ(info.offset_slack, n_tangent);
    EXPECT_EQ(info.pd_orig_offset_slack, n_tangent);
    EXPECT_LE(info.offset_slack + info.number_of_slack_variables, n_tangent + n_slack);
}

    } // namespace test
} // namespace fatrop
