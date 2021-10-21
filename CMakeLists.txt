# Set the minimum version of CMake that can be used
# To find the cmake version run
# $ cmake --version
cmake_minimum_required(VERSION 3.14)
project(FATROP)
# minimum cmake version to find numpy header
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Debug" CACHE
      STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 ${PROFILINGFLAG} -DDEBUG -DLEVEL1_DCACHE_LINE_SIZE=${L1_CACHE_LINE_SIZE} ${VERBOSEFLAG} -DPANEL_SIZE=${PANEL_SIZE} -Wall -Wextra")
set(CMAKE_C_FLAGS_DEBUG "-g -O0 ${PROFILINGFLAG} -DDEBUG -DLEVEL1_DCACHE_LINE_SIZE=${L1_CACHE_LINE_SIZE} ${VERBOSEFLAG} -DPANEL_SIZE=${PANEL_SIZE} -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DRELEASE -DLEVEL1_DCACHE_LINE_SIZE=${L1_CACHE_LINE_SIZE} -DPANEL_SIZE=${PANEL_SIZE}")
option(BUILD_DOCS "Build documentation" OFF)

find_package(blasfeo)
add_subdirectory(src)
add_subdirectory(examples)

if(BUILD_DOCS MATCHES ON)
  add_subdirectory(docs)
endif()
