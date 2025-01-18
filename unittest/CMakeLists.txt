# Option to use system Google Test (default is OFF)
option(WITH_BUILD_GTEST "Use system-installed Google Test" OFF)

if(NOT WITH_BUILD_GTEST)
  # Try to find system-installed Google Test
  find_package(GTest)
endif()

if(NOT WITH_BUILD_GTEST AND GTest_FOUND)
  message(STATUS "Using system-installed Google Test")
else()
  if(NOT WITH_BUILD_GTEST)
    message(STATUS "System Google Test not found. Falling back to fetching Google Test.")
  endif()

# Add Google Test
include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/5376968f6948923e2411081fd9372e71a59d8e77.zip
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  )
  FetchContent_MakeAvailable(googletest)
endif()

# Add Google Test
include(GoogleTest)

# Add the test executables
add_executable(vector_test linear_algebra/vector_test.cpp)
add_executable(matrix_test linear_algebra/matrix_test.cpp)
add_executable(lu_factorization_test linear_algebra/lu_factorization_test.cpp)
add_executable(jacobian_hessian_test ocp/jacobian_hessian_test.cpp)
add_executable(aug_system_solver_test ocp/aug_system_solver_test.cpp)
add_executable(problem_info_test ocp/problem_info_test.cpp)
add_executable(pd_system_solver_test ocp/pd_system_solver_test.cpp)
add_executable(ocp_test ocp/ocp_test.cpp)
add_executable(nlp_ocp_test ocp/nlp_ocp_test.cpp)

# Link the test executables with Google Test and other required libraries
target_link_libraries(vector_test gtest gtest_main blasfeo fatrop)
target_link_libraries(matrix_test gtest gtest_main blasfeo fatrop)
target_link_libraries(lu_factorization_test gtest gtest_main blasfeo fatrop)
target_link_libraries(jacobian_hessian_test gtest gtest_main blasfeo fatrop)
target_link_libraries(aug_system_solver_test gtest gtest_main blasfeo fatrop)
target_link_libraries(problem_info_test gtest gtest_main blasfeo fatrop)
target_link_libraries(pd_system_solver_test gtest gtest_main blasfeo fatrop)
target_link_libraries(ocp_test gtest gtest_main blasfeo fatrop)
target_link_libraries(nlp_ocp_test gtest gtest_main blasfeo fatrop)

# Discover and add the tests
gtest_discover_tests(vector_test)
gtest_discover_tests(matrix_test)
gtest_discover_tests(lu_factorization_test)
gtest_discover_tests(jacobian_hessian_test)
gtest_discover_tests(aug_system_solver_test)
gtest_discover_tests(problem_info_test)
gtest_discover_tests(pd_system_solver_test)
gtest_discover_tests(ocp_test)
gtest_discover_tests(nlp_ocp_test)
