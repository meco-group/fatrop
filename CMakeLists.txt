cmake_minimum_required(VERSION 3.10)
project(fatrop VERSION 1.0.0 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Generate position-independent code
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Set default build type to Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Set compilation flags for Release and Debug
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")

# Option to build tests (default is ON)
option(BUILD_TESTS "Build the tests" ON)
option(BUILD_WITH_LEGACY "Build with legacy support" ON)

# Define installation directories
set(LIBRARY_INSTALL_DIR lib)
set(ARCHIVE_INSTALL_DIR lib)
set(RUNTIME_INSTALL_DIR bin)
set(INCLUDE_INSTALL_DIR include)
set(CMAKE_INSTALL_DIR cmake)

# Add the external directory
add_subdirectory(external)

# Add BLASFEO
if(WITH_BUILD_BLASFEO)
    # BLASFEO is already available from FetchContent
    add_library(blasfeo::blasfeo ALIAS blasfeo)
else()
    find_package(blasfeo REQUIRED)
endif()


# Add include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# Add the fatrop source directory
add_subdirectory(src)

# Add the legacy directory if BUILD_WITH_LEGACY is enabled
if(BUILD_WITH_LEGACY)
    add_subdirectory(legacy)
endif()

# Add the examples directory
add_subdirectory(examples)

# Build tests if the option is enabled
if(BUILD_TESTS)
  # Enable testing
  enable_testing()

  # Add the unittest directory, but exclude it from installation
  add_subdirectory(unittest)
endif()
