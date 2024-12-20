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

# Link the test executables with Google Test and other required libraries
target_link_libraries(vector_test gtest gtest_main blasfeo)
target_link_libraries(matrix_test gtest gtest_main blasfeo)

# Discover and add the tests
gtest_discover_tests(vector_test)
gtest_discover_tests(matrix_test)
