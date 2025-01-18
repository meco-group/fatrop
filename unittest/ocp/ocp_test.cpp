#include "ocp_test_probem.hpp"
#include <fatrop/ocp/ocp_abstract.hpp>
#include <gtest/gtest.h>

using namespace fatrop::test;
class OCPTest : public ::testing::Test
{
protected:
    OcpTestProblem ocp;
    void SetUp() override
    {
        // Set up test environment
    }

    void TearDown() override
    {
        // Clean up after test
    }

    // Declare any objects you need for testing
    // FatropOCP ocp;
};

TEST_F(OCPTest, TestName)
{
    // Your test code here
    EXPECT_TRUE(true); // Replace with actual test
}

// Add more test cases as needed

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
