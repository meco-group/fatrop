#include "fatrop/common/options.hpp"
#include <gtest/gtest.h>

using namespace fatrop;
// Mock algorithm class for testing
class MockAlgo
{
public:
    void set_int_option(const Index &value) { int_option = value; }
    void set_double_option(const Scalar &value) { double_option = value; }
    void set_bool_option(const bool &value) { bool_option = value; }
    void set_string_option(const std::string &value) { string_option = value; }

    void register_options(OptionRegistry &registry)
    {
        registry.register_option("int_option", &MockAlgo::set_int_option, this);
        registry.register_option("double_option", &MockAlgo::set_double_option, this);
        registry.register_option("bool_option", &MockAlgo::set_bool_option, this);
        registry.register_option("string_option", &MockAlgo::set_string_option, this);
    }

    Index int_option = 0;
    Scalar double_option = 0.0;
    bool bool_option = false;
    std::string string_option;
};

TEST(OptionRegistryTest, RegisterAndSetOptions)
{
    OptionRegistry registry;
    MockAlgo algo;

    // Register options
    registry.register_options(algo);

    // Set and check integer option
    registry.set_option("int_option", 42);
    EXPECT_EQ(algo.int_option, 42);

    // Set and check double option
    registry.set_option("double_option", 3.14);
    EXPECT_DOUBLE_EQ(algo.double_option, 3.14);

    // Set and check bool option
    registry.set_option("bool_option", true);
    EXPECT_TRUE(algo.bool_option);

    // Set and check string option
    registry.set_option("string_option", "test");
    EXPECT_EQ(algo.string_option, "test");
}

TEST(OptionRegistryTest, SetBoolOptionWithInt)
{
    OptionRegistry registry;
    MockAlgo algo;

    registry.register_options(algo);

    // Set bool option with int
    registry.set_option("bool_option", 1);
    EXPECT_TRUE(algo.bool_option);

    registry.set_option("bool_option", 0);
    EXPECT_FALSE(algo.bool_option);
}

TEST(OptionRegistryTest, SetBoolOptionWithString)
{
    OptionRegistry registry;
    MockAlgo algo;

    registry.register_options(algo);

    // Set bool option with string
    registry.set_option("bool_option", "yes");
    EXPECT_TRUE(algo.bool_option);

    registry.set_option("bool_option", "no");
    EXPECT_FALSE(algo.bool_option);

    // Invalid string value
    EXPECT_THROW(registry.set_option("bool_option", "invalid"), std::runtime_error);
}

TEST(OptionRegistryTest, InvalidOptionType)
{
    OptionRegistry registry;
    MockAlgo algo;

    registry.register_options(algo);

    // Try to set int option with string
    EXPECT_THROW(registry.set_option("int_option", "not an int"), std::runtime_error);

    // Try to set double option with bool
    EXPECT_THROW(registry.set_option("double_option", true), std::runtime_error);

    // Try to set string option with int
    EXPECT_THROW(registry.set_option("string_option", 42), std::runtime_error);
}

TEST(OptionRegistryTest, NonExistentOption)
{
    OptionRegistry registry;
    MockAlgo algo;

    registry.register_options(algo);

    // Try to set a non-existent option
    EXPECT_THROW(registry.set_option("non_existent_option", 42), std::runtime_error);
}

TEST(OptionRegistryTest, MultipleRegistration)
{
    OptionRegistry registry;
    MockAlgo algo1, algo2;

    // Register options for both algorithms
    registry.register_options(algo1);
    registry.register_options(algo2);

    // Set option and check if it's set for both algorithms
    registry.set_option("int_option", 42);
    EXPECT_EQ(algo1.int_option, 42);
    EXPECT_EQ(algo2.int_option, 42);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}