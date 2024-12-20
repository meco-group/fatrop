cmake_minimum_required(VERSION 3.10)
project(fatrop)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Option to build tests (default is ON)
option(BUILD_TESTS "Build the tests" ON)

# Add BLASFEO
find_package(blasfeo REQUIRED)

# Add include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# Add the examples directory
add_subdirectory(examples)

# Build tests if the option is enabled
if(BUILD_TESTS)
  # Enable testing
  enable_testing()

  # Add the unittest directory
  add_subdirectory(unittest)
endif()
