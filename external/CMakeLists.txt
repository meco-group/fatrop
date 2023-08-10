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

option(BLASFEO_FROM_MAKE "Use BLASFEO from make" OFF)

if(BLASFEO_FROM_MAKE)
    # Set default values for BLASFEO header and library directories
    set(BLASFEO_HEADER_DIR "${CMAKE_CURRENT_LIST_DIR}/blasfeo/include" CACHE PATH "BLASFEO header installation directory")
    set(BLASFEO_LIB_DIR "${CMAKE_CURRENT_LIST_DIR}/blasfeo/lib" CACHE PATH "BLASFEO library directory")

    # Create an imported static target for blasfeo
    add_library(blasfeo STATIC IMPORTED GLOBAL)

    # Set the include directories and library file for blasfeo
    set_target_properties(blasfeo PROPERTIES
        IMPORTED_LOCATION ${BLASFEO_LIB_DIR}/libblasfeo.a
        INTERFACE_INCLUDE_DIRECTORIES ${BLASFEO_HEADER_DIR}
    )

    # Print message
    message(STATUS "Using BLASFEO from make. Make sure you compiled BLASFEO manually.")
    message(STATUS "BLASFEO Header Directory: ${BLASFEO_HEADER_DIR}")
    message(STATUS "BLASFEO Library Directory: ${BLASFEO_LIB_DIR}")
else()
    message(STATUS "using BLASFEO with cmake from git submodule.")
    set(BLASFEO_HEADERS_INSTALLATION_DIRECTORY "include/blasfeo/include" CACHE STRING "")
    set(TARGET ${BLASFEO_TARGET} CACHE STRING "Set CPU architecture target" FORCE)
    set(MF "PANELMAJ" CACHE STRING "Matrix format" FORCE)
    set(BLAS_API OFF CACHE BOOL "Compile BLAS API" FORCE)
    set(BLASFEO_REF_API ON CACHE BOOL "REF")
    set(BLASFEO_EXAMPLES OFF CACHE BOOL "Build blasfeo examples")
    add_subdirectory(blasfeo)
endif()

if(WIN32)
    add_subdirectory(dlfcn-win32)
endif(WIN32)