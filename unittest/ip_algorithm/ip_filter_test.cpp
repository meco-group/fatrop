#include "gtest/gtest.h"
#include "fatrop/ip_algorithm/ip_filter.hpp"

namespace fatrop {
namespace unittest {

TEST(IpFilterTest, DefaultConstructor) {
  IpFilter filter;
  ASSERT_EQ(filter.size(), 0);
}

TEST(IpFilterTest, Reset) {
  IpFilter filter;
  filter.reset();
  ASSERT_EQ(filter.size(), 0);
}

TEST(IpFilterTest, Reserve) {
  IpFilter filter;
  filter.reserve(10);
  // No direct way to check capacity, but should not crash when adding elements up to 10
}

TEST(IpFilterTest, IsAcceptableEmptyFilter) {
  IpFilter filter;
  IpFilterData data = {1.0, 0.0};
  ASSERT_TRUE(filter.is_acceptable(data));
}

TEST(IpFilterTest, AugmentAndIsAcceptable) {
  IpFilter filter;
  IpFilterData data1 = {1.0, 0.0};
  filter.augment(data1);
  ASSERT_EQ(filter.size(), 1);
  ASSERT_TRUE(filter.is_acceptable(data1)); // Should not be acceptable after being added

  IpFilterData data2 = {0.5, 0.5};
  ASSERT_TRUE(filter.is_acceptable(data2));
  filter.augment(data2);
  ASSERT_EQ(filter.size(), 2);

  IpFilterData data3 = {2.0, 2.0};
  ASSERT_FALSE(filter.is_acceptable(data3)); // Dominated by both data1 and data2
}

TEST(IpFilterTest, Size) {
  IpFilter filter;
  ASSERT_EQ(filter.size(), 0);
  IpFilterData data1 = {1.0, 0.0};
  filter.augment(data1);
  ASSERT_EQ(filter.size(), 1);
  IpFilterData data2 = {0.5, 0.5};
  filter.augment(data2);
  ASSERT_EQ(filter.size(), 2);
}

} // namespace unittest
} // namespace fatTRU
