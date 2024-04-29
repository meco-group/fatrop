#[[
Fatrop - A fast trajectory optimization solver
Copyright (C) 2022, 2023 Lander Vanroye, KU Leuven. All rights reserved.

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

option(BUILD_WITH_BLASFEO "Build BLASFEO with cmake from git submodule" ON)

if (BUILD_WITH_BLASFEO)
    message(STATUS "using BLASFEO with cmake from git submodule.")
    set(BLASFEO_HEADERS_INSTALLATION_DIRECTORY "include/blasfeo/include" CACHE STRING "")
    set(TARGET ${BLASFEO_TARGET} CACHE STRING "Set CPU architecture target" FORCE)
    set(MF "PANELMAJ" CACHE STRING "Matrix format" FORCE)
    set(BLAS_API OFF CACHE BOOL "Compile BLAS API" FORCE)
    set(BLASFEO_REF_API ON CACHE BOOL "REF")
    set(BLASFEO_EXAMPLES OFF CACHE BOOL "Build blasfeo examples")
    add_subdirectory(blasfeo)

else()
    set(BLASFEO_INSTALL_DIR "/opt/blasfeo" CACHE PATH "Path to blasfeo installation directory")
    # find_library(blasfeo NAMES blasfeo PATHS ${BLASFEO_INSTALL_DIR}/lib)

    find_path(BLASFEO_INCLUDE_DIR NAMES "blasfeo.h" PATHS ${BLASFEO_INSTALL_DIR}/include)
    find_library(BLASFEO_LIBRARY NAMES blasfeo HINTS ${BLASFEO_INSTALL_DIR}/lib)

    # Create an imported static target for blasfeo
    add_library(blasfeo STATIC IMPORTED GLOBAL)
    # Set the path to the library file
    set_target_properties(blasfeo PROPERTIES IMPORTED_LOCATION ${BLASFEO_LIBRARY})

    # Set the path to the include directory
    set_target_properties(blasfeo PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${BLASFEO_INCLUDE_DIR})
    
endif()

if(WITH_PYTHON)
    # add pybind11 for python bindings
    add_subdirectory(pybind11)
endif()

if(WIN32)
    add_subdirectory(dlfcn-win32)
endif(WIN32)