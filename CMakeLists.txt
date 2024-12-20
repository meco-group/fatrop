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

# Define source files for the fatrop library
set(FATROP_SOURCES
    src/linear_algebra/lu_factorization.cpp
    # Add other source files here
)

# Create the fatrop shared library
add_library(fatrop SHARED ${FATROP_SOURCES})

# Link fatrop with blasfeo
target_link_libraries(fatrop PUBLIC blasfeo)

# Add the examples directory
add_subdirectory(examples)

# Build tests if the option is enabled
if(BUILD_TESTS)
  # Enable testing
  enable_testing()

  # Add the unittest directory
  add_subdirectory(unittest)
endif()
