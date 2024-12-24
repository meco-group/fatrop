#include "fatrop/common/vector_index.hpp"
#include <gtest/gtest.h>
#include <vector>

using namespace fatrop;

// Test fixture for VecIndex and VecIndexBlock
template <typename T>
class VecIndexTest : public ::testing::Test {
protected:
    const std::vector<T> test_data = {1, 2, 3, 4, 5};
    VecIndexView<T> vec_index_view;

    VecIndexTest() : vec_index_view(test_data) {}
};

using TestedTypes = ::testing::Types<int, long, size_t>;
TYPED_TEST_SUITE(VecIndexTest, TestedTypes);

TYPED_TEST(VecIndexTest, OperatorBrackets) {
    EXPECT_EQ(this->vec_index_view[2], 3);
    EXPECT_EQ(this->vec_index_view[0], 1);
    EXPECT_EQ(this->vec_index_view[4], 5);
}

TYPED_TEST(VecIndexTest, Size) {
    EXPECT_EQ(this->vec_index_view.m(), 5);
}

TYPED_TEST(VecIndexTest, Block) {
    auto block = this->vec_index_view.block(1, 3);
    EXPECT_EQ(block.m(), 3);
    EXPECT_EQ(block[0], 2);
    EXPECT_EQ(block[1], 3);
    EXPECT_EQ(block[2], 4);
}

TYPED_TEST(VecIndexTest, Sum) {
    EXPECT_EQ(sum(this->vec_index_view), 15);
}

TYPED_TEST(VecIndexTest, VecIndexBlockTest) {
    auto block = this->vec_index_view.block(1, 3);
    EXPECT_EQ(block.m(), 3);
    EXPECT_EQ(block[0], 2);
    EXPECT_EQ(block[1], 3);
    EXPECT_EQ(block[2], 4);

    // Test nested blocks
    auto nested_block = block.block(1, 2);
    EXPECT_EQ(nested_block.m(), 2);
    EXPECT_EQ(nested_block[0], 3);
    EXPECT_EQ(nested_block[1], 4);
}

// Test fixture for VecIndexView
template <typename T>
class IndexVecViewTest : public ::testing::Test {
protected:
    std::vector<T> test_data = {1, 2, 3, 4, 5};
    VecIndexView<T> index_vec_view;

    IndexVecViewTest() : index_vec_view(test_data, 1, 3) {}
};

TYPED_TEST_SUITE(IndexVecViewTest, TestedTypes);

TYPED_TEST(IndexVecViewTest, Constructor) {
    EXPECT_EQ(this->index_vec_view.m(), 3);
    EXPECT_EQ(this->index_vec_view[0], 2);
    EXPECT_EQ(this->index_vec_view[1], 3);
    EXPECT_EQ(this->index_vec_view[2], 4);
}

TYPED_TEST(IndexVecViewTest, Block) {
    auto block = this->index_vec_view.block(1, 2);
    EXPECT_EQ(block.m(), 2);
    EXPECT_EQ(block[0], 3);
    EXPECT_EQ(block[1], 4);
}

TYPED_TEST(VecIndexTest, FullVectorConstructor) {
    VecIndexView<TypeParam> full_view(this->test_data);
    EXPECT_EQ(full_view.m(), 5);
    for (size_t i = 0; i < 5; ++i) {
        EXPECT_EQ(full_view[i], this->test_data[i]);
    }
}

TYPED_TEST(VecIndexTest, ConversionToStdVector) {
    std::vector<TypeParam> converted = this->vec_index_view;
    EXPECT_EQ(converted.size(), 5);
    for (size_t i = 0; i < 5; ++i) {
        EXPECT_EQ(converted[i], this->test_data[i]);
    }
}

TYPED_TEST(VecIndexTest, EmptyVector) {
    std::vector<TypeParam> empty_data;
    VecIndexView<TypeParam> empty_view(empty_data);
    EXPECT_EQ(empty_view.m(), 0);
}

// Note: This test assumes that the VecIndexView class doesn't perform bounds checking.
// If it does, you may need to modify this test or remove it.
TYPED_TEST(VecIndexTest, OutOfBoundsAccess) {
    EXPECT_NO_THROW(this->vec_index_view[4]); // Last valid index
    EXPECT_NO_THROW(this->vec_index_view[5]); // Out of bounds
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
