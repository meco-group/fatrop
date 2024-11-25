#include <gtest/gtest.h>
#include <stdexcept>
#include "fatrop/solver/FatropOptions.hpp"  // Replace with the actual file name containing your code

namespace fatrop {

// Test OptionBase and derived classes functionality
class OptionTests : public ::testing::Test {
protected:
    FatropOptions options;  // Create an instance of FatropOptions to test
    FatropOptionsRegistry registry{options};  // Create the registry with the options
};

// Test valid setting for IntOption
TEST_F(OptionTests, TestSetIntOption) {
    OptionValueVariant new_value(200);  // An int value to set
    registry.set("max_iter", new_value);

    EXPECT_EQ(options.max_iter.get(), 200);  // Validate the value was set correctly
}

// Test valid setting for DoubleOption
TEST_F(OptionTests, TestSetDoubleOption) {
    OptionValueVariant new_value(1e-5);  // A double value to set
    registry.set("tol", new_value);

    EXPECT_EQ(options.tol.get(), 1e-5);  // Validate the value was set correctly
}

// Test valid setting for BoolOption with "yes"
TEST_F(OptionTests, TestSetBoolOptionYes) {
    OptionValueVariant new_value("yes");
    registry.set("recalc_y", new_value);

    EXPECT_TRUE(options.recalc_y.get());  // Validate the boolean value was set to true
}

// Test valid setting for BoolOption with "no"
TEST_F(OptionTests, TestSetBoolOptionNo) {
    OptionValueVariant new_value("no");
    registry.set("recalc_y", new_value);

    EXPECT_FALSE(options.recalc_y.get());  // Validate the boolean value was set to false
}

// Test setting out-of-bounds value for IntOption (lower bound)
TEST_F(OptionTests, TestSetIntOptionOutOfBoundsLower) {
    OptionValueVariant new_value(-1);  // Invalid value (out of bounds for max_iter)

    EXPECT_THROW(registry.set("max_iter", new_value), std::invalid_argument);  // Should throw an exception
}

// Test setting out-of-bounds value for DoubleOption (lower bound)
TEST_F(OptionTests, TestSetDoubleOptionOutOfBoundsLower) {
    OptionValueVariant new_value(-1e-7);  // Invalid value (out of bounds for tol)

    EXPECT_THROW(registry.set("tol", new_value), std::invalid_argument);  // Should throw an exception
}

// Test setting invalid type (wrong type, e.g. string for IntOption)
TEST_F(OptionTests, TestSetInvalidTypeForIntOption) {
    OptionValueVariant new_value("wrong type");  // Invalid type (string instead of int)

    EXPECT_THROW(registry.set("max_iter", new_value), std::invalid_argument);  // Should throw an exception
}

// Test setting invalid string value for BoolOption
TEST_F(OptionTests, TestSetInvalidStringForBoolOption) {
    OptionValueVariant new_value("maybe");  // Invalid string for a boolean option

    EXPECT_THROW(registry.set("recalc_y", new_value), std::invalid_argument);  // Should throw an exception
}

// Test setting value for an option that doesn't exist
TEST_F(OptionTests, TestSetNonExistentOption) {
    OptionValueVariant new_value(42);  // Arbitrary value

    EXPECT_THROW(registry.set("non_existent_option", new_value), std::invalid_argument);  // Should throw an exception
}

}  // namespace fatrop

// Main function to run the tests
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
