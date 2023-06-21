#[[
Fatrop - A fast trajectory optimization solver
Copyright (C) 2022, 2023 Lander Vanroye <lander.vanroye@kuleuven.be>

This file is part of Fatrop.

Fatrop is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Fatrop is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with Fatrop.  If not, see <http://www.gnu.org/licenses/>.]]

# Set the minimum version of CMake that can be used
# To find the cmake version run
# $ cmake --version
cmake_minimum_required(VERSION 3.14)
set(VERSION 0.0.1)
project(FATROP VERSION "${VERSION}" LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

# minimum cmake version to find numpy header
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Release" CACHE
      STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

if (WIN32)
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
  set(CMAKE_CXX_FLAGS_DEBUG "-O0 ${PROFILINGFLAG} -DDEBUG -DLEVEL1_DCACHE_LINE_SIZE=${L1_CACHE_LINE_SIZE} ${VERBOSEFLAG} -DPANEL_SIZE=${PANEL_SIZE} -Wall")
  set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DRELEASE -DLEVEL1_DCACHE_LINE_SIZE=${L1_CACHE_LINE_SIZE} -DPANEL_SIZE=${PANEL_SIZE}")
else()
  set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 ${PROFILINGFLAG} -DDEBUG -DLEVEL1_DCACHE_LINE_SIZE=${L1_CACHE_LINE_SIZE} ${VERBOSEFLAG} -DPANEL_SIZE=${PANEL_SIZE} -Wall -Wextra -Wno-unused-parameter")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DRELEASE -DLEVEL1_DCACHE_LINE_SIZE=${L1_CACHE_LINE_SIZE} -DPANEL_SIZE=${PANEL_SIZE}")
endif(WIN32)

option(BUILD_DOCS "Build documentation" OFF)

# TODO write this is a more modern way
set(ENABLE_MULTITHREADING OFF CACHE BOOL "Enable multithreading")

if(${ENABLE_MULTITHREADING})
  set(CMAKE_CXX_STANDARD 20)
  set(CMAKE_CXX_COMPILER "g++-11")
  add_definitions(-DENABLE_MULTITHREADING)
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=gnu++20 -fopenmp")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=gnu++20 -fopenmp")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread")
endif()
set(BLASFEO_TARGET "X64_AUTOMATIC" CACHE STRING "BLASFEO Target architecture")
set(BLASFEO_SRC_DIR ${PROJECT_SOURCE_DIR}/external/blasfeo)
message(STATUS "Target: BLASFEO is ${BLASFEO_TARGET}")
add_subdirectory(external)
add_subdirectory(fatrop)

if(BUILD_DOCS MATCHES ON)
  add_subdirectory(docs)
endif()