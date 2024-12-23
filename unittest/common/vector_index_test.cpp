#include "fatrop/common/vector_index.hpp"
#include <gtest/gtest.h>
#include <vector>

using namespace fatrop;

// Test fixture for VecIndex and VecIndexBlock
template <typename T>
class VecIndexTest : public ::testing::Test {
protected:
    std::vector<T> test_data = {1, 2, 3, 4, 5};
    VecIndexAllocated<T> vec_index_allocated;

    VecIndexTest() : vec_index_allocated(test_data) {}
};

using TestedTypes = ::testing::Types<int, long, size_t>;
TYPED_TEST_SUITE(VecIndexTest, TestedTypes);

TYPED_TEST(VecIndexTest, OperatorBrackets) {
    EXPECT_EQ(this->vec_index_allocated[2], 3);
    EXPECT_EQ(this->vec_index_allocated[0], 1);
    EXPECT_EQ(this->vec_index_allocated[4], 5);
}

TYPED_TEST(VecIndexTest, Size) {
    EXPECT_EQ(this->vec_index_allocated.m(), 5);
}

TYPED_TEST(VecIndexTest, Block) {
    auto block = this->vec_index_allocated.block(1, 3);
    EXPECT_EQ(block.m(), 3);
    EXPECT_EQ(block[0], 2);
    EXPECT_EQ(block[1], 3);
    EXPECT_EQ(block[2], 4);
}

TYPED_TEST(VecIndexTest, Sum) {
    EXPECT_EQ(sum(this->vec_index_allocated), 15);
}

TYPED_TEST(VecIndexTest, VecIndexBlockTest) {
    auto block = this->vec_index_allocated.block(1, 3);
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

// Test fixture for IndexVecView
template <typename T>
class IndexVecViewTest : public ::testing::Test {
protected:
    std::vector<T> test_data = {1, 2, 3, 4, 5};
    VecIndexAllocated<T> vec_index_allocated;
    IndexVecView<T> index_vec_view;

    IndexVecViewTest() : vec_index_allocated(test_data), index_vec_view(vec_index_allocated, 1, 3) {}
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

TYPED_TEST(IndexVecViewTest, AssignmentOperator) {
    std::vector<TypeParam> new_data = {10, 20, 30};
    VecIndexAllocated<TypeParam> new_vec(new_data);
    this->index_vec_view = new_vec;
    EXPECT_EQ(this->index_vec_view[0], 10);
    EXPECT_EQ(this->index_vec_view[1], 20);
    EXPECT_EQ(this->index_vec_view[2], 30);
}

// Test fixture for VecIndexAllocated
template <typename T>
class VecIndexAllocatedTest : public ::testing::Test {
protected:
    std::vector<T> test_data = {1, 2, 3, 4, 5};
};

TYPED_TEST_SUITE(VecIndexAllocatedTest, TestedTypes);

TYPED_TEST(VecIndexAllocatedTest, ConstructorRValue) {
    VecIndexAllocated<TypeParam> vec(std::vector<TypeParam>{1, 2, 3, 4, 5});
    EXPECT_EQ(vec.m(), 5);
    EXPECT_EQ(vec[2], 3);
}

TYPED_TEST(VecIndexAllocatedTest, ConstructorConstRef) {
    VecIndexAllocated<TypeParam> vec(this->test_data);
    EXPECT_EQ(vec.m(), 5);
    EXPECT_EQ(vec[2], 3);
}

// TYPED_TEST(VecIndexAllocatedTest, EmptyConstructor) {
//     VecIndexAllocated<TypeParam> vec;
//     EXPECT_EQ(vec.m(), 0);
// }

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
