# Set the minimum version of CMake that can be used
# To find the cmake version run
# $ cmake --version
cmake_minimum_required(VERSION 3.14)
set(VERSION 0.0.1)
project(FATROP VERSION "${VERSION}" LANGUAGES C CXX Fortran)

set(CMAKE_CXX_STANDARD 11)

# include(FetchContent)
# FetchContent_Declare(
#   googletest
#   URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
# )
# For Windows: Prevent overriding the parent project's compiler/linker settings
# set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
# FetchContent_MakeAvailable(googletest)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

# minimum cmake version to find numpy header
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Debug" CACHE
      STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()
execute_process(COMMAND "getconf" "LEVEL1_DCACHE_LINESIZE" OUTPUT_VARIABLE L1_CACHE_LINE_SIZE OUTPUT_STRIP_TRAILING_WHITESPACE)
set(CACHE_LINE_SIZE "${L1_CACHE_LINE_SIZE}" CACHE STRING "cache line size used for cache aligning")
math(EXPR L1_CACHE_LINE_SIZE "${CACHE_LINE_SIZE}")

# set(CMAKE_CXX_FLAGS_DEBUG "-pg -fprofile-arcs -ftest-coverage -g -O0 ${PROFILINGFLAG} -DDEBUG -DLEVEL1_DCACHE_LINE_SIZE=${L1_CACHE_LINE_SIZE} ${VERBOSEFLAG} -DPANEL_SIZE=${PANEL_SIZE} -Wall -Wextra -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 ${PROFILINGFLAG} -DDEBUG -DLEVEL1_DCACHE_LINE_SIZE=${L1_CACHE_LINE_SIZE} ${VERBOSEFLAG} -DPANEL_SIZE=${PANEL_SIZE} -Wall -Wextra -Wno-unused-parameter")
set(CMAKE_C_FLAGS_DEBUG "-g -O0 ${PROFILINGFLAG} -DDEBUG -DLEVEL1_DCACHE_LINE_SIZE=${L1_CACHE_LINE_SIZE} ${VERBOSEFLAG} -DPANEL_SIZE=${PANEL_SIZE} -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DRELEASE -DLEVEL1_DCACHE_LINE_SIZE=${L1_CACHE_LINE_SIZE} -DPANEL_SIZE=${PANEL_SIZE}")
# set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fno-split-loops -fno-loop-unroll-and-jam -fno-tree-loop-vectorize -fno-split-loops -fno-align-loops -fno-loop-interchange -march=native -DRELEASE -DLEVEL1_DCACHE_LINE_SIZE=${L1_CACHE_LINE_SIZE} -DPANEL_SIZE=${PANEL_SIZE}")
option(BUILD_DOCS "Build documentation" OFF)

# TODO write this is a more modern way
set(ENABLE_MULTITHREADING OFF CACHE BOOL "Enable multithreading")
if(${ENABLE_MULTITHREADING})
  set(CMAKE_CXX_STANDARD 20)
  set(CMAKE_C_COMPILER "gcc-11")
  set(CMAKE_CXX_COMPILER "g++-11")
  add_definitions(-DENABLE_MULTITHREADING)
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=gnu++20")
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -std=gnu++20")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=gnu++20")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread")
endif()

# include(FindMUMPS.cmake)
# find_package(MUMPS)
find_package(blasfeo)
add_subdirectory(fatrop)
# add_subdirectory(examples)
# add_subdirectory(SWIG)
# add_subdirectory(test)

if(BUILD_DOCS MATCHES ON)
  add_subdirectory(docs)
endif()
