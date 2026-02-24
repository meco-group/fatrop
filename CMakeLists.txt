cmake_minimum_required(VERSION 3.14)
project(fatrop VERSION 1.0.0 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Generate position-independent code
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Set default build type to Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release")
endif()

# Set compilation flags for Release and Debug
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")

# Option to build tests (default is ON)
option(BUILD_TESTS "Build the tests" ON)
option(BUILD_EXAMPLES "Build the examples" ON)
option(BUILD_WITH_C_INTERFACE "Build with C interface support" ON)

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

# Add the c_interface directory if BUILD_WITH_C_INTERFACE is enabled
if(BUILD_WITH_C_INTERFACE)
    add_subdirectory(c_interface)
endif()

# Add the examples directory
if(BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

# Build tests if the option is enabled
if(BUILD_TESTS)
  # Enable testing
  enable_testing()

  # Add the unittest directory, but exclude it from installation
  add_subdirectory(unittest)
endif()
