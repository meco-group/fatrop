#include "fatrop/common/vector_index.hpp"
#include <gtest/gtest.h>
#include <vector>

using namespace fatrop;

class IndexVectorTest : public ::testing::Test
{
protected:
    VecIndexAllocated<int> test_data = std::vector<int>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
};
TEST_F(IndexVectorTest, SumTest)
{
    std::cout << sum(test_data) << std::endl;
    std::cout << sum(test_data.block(1, 8)) << std::endl;
    // std::vector<int> test2(10);
    // test2 = test_data;
    // for(int i = 0; i < 10; i++) {
    //     std::cout << test2[i] << std::endl;
    // }

    // EXPECT_EQ(sum(test_data), 45);
}

// TEST_F(IndexVectorTest, SumTest) {
//     EXPECT_EQ(sum(test_data), 45);
// }

// TEST_F(IndexVectorTest, SumBlockTest) {
//     EXPECT_EQ(sum(block(test_data, 1, 8)), 36);
// }

// TEST_F(IndexVectorTest, AssignTest) {
//     EXPECT_EQ(sum(block(test_data, 1, 8)), 36);
// }

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
