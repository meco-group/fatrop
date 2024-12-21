#include "fatrop/ocp/jacobian.hpp"
#include "fatrop/ocp/hessian.hpp"
#include "fatrop/ocp/dims.hpp"
#include "fatrop/ocp/type.hpp"
#include "fatrop/context/context.hpp"
#include <gtest/gtest.h>
#include <vector>

using namespace fatrop;

TEST(JacobianTest, ConstructorTest) {
    // Create OcpDims object
    int K = 5;  // Number of stages
    std::vector<Index> nx = {2, 2, 2, 2, 2, 2};  // State dimensions for each stage
    std::vector<Index> nu = {1, 1, 1, 1, 1};     // Input dimensions for each stage
    std::vector<Index> ng = {1, 0, 0, 0, 0, 2};  // Equality constraints for each stage
    std::vector<Index> ng_ineq = {1, 0, 0, 2, 0, 0};  // Inequality constraints for each stage

    OcpDims dims(K, nx, nu, ng, ng_ineq);

    // Check if Jacobian object can be constructed without throwing an exception
    EXPECT_NO_THROW({
        Jacobian<OcpType> jacobian(dims);
    });
}

TEST(HessianTest, ConstructorTest) {
    // Create OcpDims object
    int K = 5;  // Number of stages
    std::vector<Index> nx = {2, 2, 2, 2, 2, 2};  // State dimensions for each stage
    std::vector<Index> nu = {1, 1, 1, 1, 1};     // Input dimensions for each stage
    std::vector<Index> ng = {1, 0, 0, 0, 0, 2};  // Equality constraints for each stage
    std::vector<Index> ng_ineq = {1, 0, 0, 2, 0, 0};  // Inequality constraints for each stage

    OcpDims dims(K, nx, nu, ng, ng_ineq);

    // Check if Hessian object can be constructed without throwing an exception
    EXPECT_NO_THROW({
        Hessian<OcpType> hessian(dims);
    });
}

TEST(JacobianTest, AssertionViolationTest) {
    // Create OcpDims object with invalid dimensions
    int K = 3;  // Number of stages
    std::vector<Index> nx = {2, 2, 2, 2};  // State dimensions for each stage
    std::vector<Index> nu = {1, 1, 1};     // Input dimensions for each stage
    std::vector<Index> ng = {4, 3, 3, 2};  // Equality constraints for each stage (violates assertion)
    std::vector<Index> ng_ineq = {0, 0, 0, 0};  // Inequality constraints for each stage

    OcpDims dims(K, nx, nu, ng, ng_ineq);

    // Check if Jacobian constructor throws an exception due to assertion violation
    EXPECT_ANY_THROW({
        Jacobian<OcpType> jacobian(dims);
    });  // Assuming fatrop_assert throws a std::runtime_error
}
