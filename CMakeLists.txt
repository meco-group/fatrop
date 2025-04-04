#[[
Fatrop - A fast trajectory optimization solver
 Copyright (C) 2022 - 2024 Lander Vanroye, KU Leuven. All rights reserved.

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

if(WIN32)
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif(WIN32)

set(MARCH_NATIVE ON CACHE BOOL "Compile with -march-native (if applicable)")

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  set(CMAKE_CXX_FLAGS_DEBUG "-DEBUG -O0 ${PROFILINGFLAG} -DDEBUG -DLEVEL1_DCACHE_LINE_SIZE=${L1_CACHE_LINE_SIZE} ${VERBOSEFLAG} -DPANEL_SIZE=${PANEL_SIZE} -Wall")
  set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DRELEASE -DLEVEL1_DCACHE_LINE_SIZE=${L1_CACHE_LINE_SIZE} -DPANEL_SIZE=${PANEL_SIZE}")
else()
  set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 ${PROFILINGFLAG} -DDEBUG -DLEVEL1_DCACHE_LINE_SIZE=${L1_CACHE_LINE_SIZE} ${VERBOSEFLAG} -DPANEL_SIZE=${PANEL_SIZE} -Wall -Wextra -Wno-unused-parameter")
  set(MARCH_FLAG "")
  if(${MARCH_NATIVE})
    set(MARCH_FLAG "-march=native ")
  endif()
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 ${MARCH_FLAG}-DRELEASE -DLEVEL1_DCACHE_LINE_SIZE=${L1_CACHE_LINE_SIZE} -DPANEL_SIZE=${PANEL_SIZE}")
endif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")

# Specify correct version of C Runtime Libary for MSVC
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
endif()

option(BUILD_DOCS "Build documentation" OFF)
option(BUILD_EXECUTABLES "Build executables" ON)

# TODO write this is a more modern way
set(ENABLE_MULTITHREADING OFF CACHE BOOL "Enable multithreading")

if(${ENABLE_MULTITHREADING})
  if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    message(SEND_ERROR "Multithreading not yet supported when using MSVC")
  endif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")

  set(CMAKE_CXX_STANDARD 20)

  # set(CMAKE_CXX_COMPILER "g++-11")
  add_definitions(-DENABLE_MULTITHREADING)
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=gnu++2a -fopenmp")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=gnu++2a -fopenmp")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread")
endif()

set(BLASFEO_TARGET "X64_AUTOMATIC" CACHE STRING "BLASFEO Target architecture")
set(BLASFEO_SRC_DIR ${PROJECT_SOURCE_DIR}/external/blasfeo)
message(STATUS "Target: BLASFEO is ${BLASFEO_TARGET}")
add_subdirectory(external)

if(${BLASFEO_REF_API})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBLASFEO_REF_API")
endif()

add_subdirectory(fatrop)

if(BUILD_DOCS MATCHES ON)
  add_subdirectory(docs)
endif()

set(WITH_SPECTOOL OFF CACHE BOOL "Build problem specification tools, built on top of casadi")
SET(INSTALL_TARGETS fatrop)

if(${WITH_SPECTOOL})
  # casadi is required for the specification tools, find using findCasadi.cmake
  include(${CMAKE_CURRENT_LIST_DIR}/cmake/FindCasadi.cmake)
  # check if casadi/core/function_internal.hpp is available in the casadi path (CASADI_INCLUDE_DIR is set by findCasadi.cmake, and can contain multiple paths)
  message (STATUS "CASADI_INCLUDE_DIR is ${CASADI_INCLUDE_DIR}")
  if(NOT EXISTS "${CASADI_INCLUDE_DIR}/casadi/core/function_internal.hpp")
    message(SEND_ERROR "${CASADI_INCLUDE_DIR}/casadi/core/function_internal.hpp not found, please install casadi with cmake option -DINSTALL_INTERNAL_HEADERS=ON")
  endif()

  set(INSTALL_TARGETS ${INSTALL_TARGETS} spectool)
  # add the WITH_SPECTOOL compiler flag
  add_definitions(-DWITH_SPECTOOL)
  add_subdirectory(fatrop/spectool)
endif()

# set(WITH_FATROP_CASADI ON CACHE BOOL "Build problem specification tools, built on top of casadi")
# if(${WITH_FATROP_CASADI})
# # casadi is required for the specification tools, find using findCasadi.cmake
# include(${CMAKE_CURRENT_LIST_DIR}/cmake/FindCasadi.cmake)
# add_subdirectory(fatrop_casadi)
# endif()
set(WITH_PYTHON OFF CACHE BOOL "Build python bindings (fatropy)")

if(WITH_PYTHON)
  if(CMAKE_VERSION VERSION_LESS "3.18")
    find_package(Python REQUIRED COMPONENTS Interpreter Development)
  else()
    find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
  endif()
  find_package(pybind11 CONFIG)
  if(NOT pybind11_FOUND)
    add_subdirectory(external/pybind11)
  endif()
  add_subdirectory(fatropy)
endif()

if(BUILD_WITH_BLASFEO)
  set(INSTALL_TARGETS ${INSTALL_TARGETS} blasfeo)
endif()

INSTALL(TARGETS ${INSTALL_TARGETS}
  EXPORT fatropConfig
  LIBRARY DESTINATION ${LIBRARY_INSTALL_DIR}
  ARCHIVE DESTINATION ${LIBRARY_INSTALL_DIR}
  INCLUDES DESTINATION ${INCLUDE_INSTALL_DIR}
  RUNTIME DESTINATION ${RUNTIME_INSTALL_DIR}
)
INSTALL(EXPORT fatropConfig NAMESPACE fatrop:: DESTINATION cmake)

# create an option for compiling the tests
option(BUILD_TESTS "Build tests" OFF)
# if the option is set to ON, add the tests directory
if(BUILD_TESTS)
  enable_testing()
  add_subdirectory(unittest)
endif()


