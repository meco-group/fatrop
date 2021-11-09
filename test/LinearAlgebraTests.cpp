#include <iostream>
#include <memory>
#include "Fatrop.hpp"
#include "FatropLinearAlgebraEigen.hpp"
#include "FatropDebugTools.hpp"
#include <gtest/gtest.h>
using namespace fatrop;
using namespace std;
TEST(HelloTest, BasicAssertions)
{
    // Expect two strings not to be equal.
    EXPECT_STRNE("hello", "world");
    // Expect equality.
    EXPECT_EQ(7 * 6, 42);
}
class test_container : public fatrop_memory_allocator
{
};
int main(int argc, char** argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}