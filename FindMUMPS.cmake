# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

#[=======================================================================[.rst:
FindMUMPS
---------

Finds the MUMPS library.
Note that MUMPS generally requires SCALAPACK and LAPACK as well.
PORD is always used, in addition to the optional Scotch + METIS.

COMPONENTS
  s d c z   list one or more. Default is "s d"
  mpiseq    for -Dparallel=false, a stub MPI & Scalapack library

Result Variables
^^^^^^^^^^^^^^^^

MUMPS_LIBRARIES
  libraries to be linked

MUMPS_INCLUDE_DIRS
  dirs to be included

MUMPS_HAVE_OPENMP
  MUMPS is using OpenMP and thus user programs must link OpenMP as well.

MUMPS_HAVE_Scotch
  MUMPS is using Scotch/METIS and thus user programs must link Scotch/METIS as well.

#]=======================================================================]

set(MUMPS_LIBRARY)  # don't endlessly append
set(CMAKE_REQUIRED_FLAGS)

include(CheckSourceCompiles)

# --- functions

function(mumps_check)

get_property(enabled_langs GLOBAL PROPERTY ENABLED_LANGUAGES)
if(NOT Fortran IN_LIST enabled_langs)
  set(MUMPS_links true)
  return()
endif()

if(NOT mpiseq IN_LIST MUMPS_FIND_COMPONENTS)
  find_package(MPI COMPONENTS C Fortran)
  find_package(SCALAPACK)
endif()

find_package(LAPACK)

set(CMAKE_REQUIRED_INCLUDES ${MUMPS_INCLUDE_DIR} ${SCALAPACK_INCLUDE_DIRS} ${LAPACK_INCLUDE_DIRS} ${MPI_Fortran_INCLUDE_DIRS} ${MPI_C_INCLUDE_DIRS})
set(CMAKE_REQUIRED_LIBRARIES ${MUMPS_LIBRARY} ${SCALAPACK_LIBRARIES} ${LAPACK_LIBRARIES} ${MPI_Fortran_LIBRARIES} ${MPI_C_LIBRARIES})

# MUMPS doesn't set any distinct symbols or procedures if OpenMP was linked,
# so we do this indirect test to see if MUMPS needs OpenMP to link.
set(omp_src "program test_omp
implicit none (type, external)
external :: mumps_ana_omp_return, MUMPS_ICOPY_32TO64_64C
call mumps_ana_omp_return()
call MUMPS_ICOPY_32TO64_64C()
end program")

check_source_compiles(Fortran ${omp_src} mumps_openmp_test)

if(NOT mumps_openmp_test)
  find_package(OpenMP COMPONENTS C Fortran)
  list(APPEND CMAKE_REQUIRED_FLAGS ${OpenMP_Fortran_FLAGS} ${OpenMP_C_FLAGS})
  list(APPEND CMAKE_REQUIRED_INCLUDES ${OpenMP_Fortran_INCLUDE_DIRS} ${OpenMP_C_INCLUDE_DIRS})
  list(APPEND CMAKE_REQUIRED_LIBRARIES ${OpenMP_Fortran_LIBRARIES} ${OpenMP_C_LIBRARIES})

  check_source_compiles(Fortran ${omp_src} MUMPS_HAVE_OPENMP)
endif()

# check if Scotch linked
find_package(Scotch COMPONENTS parallel ESMUMPS)
# METIS is required when using Scotch
if(Scotch_FOUND)
  find_package(METIS COMPONENTS parallel)
endif()

if(METIS_FOUND AND Scotch_FOUND)
  list(APPEND CMAKE_REQUIRED_INCLUDES ${Scotch_INCLUDE_DIRS} ${METIS_INCLUDE_DIRS})
  list(APPEND CMAKE_REQUIRED_LIBRARIES ${Scotch_LIBRARIES} ${METIS_LIBRARIES})

  check_source_compiles(Fortran
  "program test_scotch
  implicit none (type, external)
  external :: mumps_dgraphinit
  call mumps_dgraphinit()
  end program"
  MUMPS_HAVE_Scotch
  )
endif()


foreach(i s d)

check_source_compiles(Fortran
  "program test_mumps
  implicit none (type, external)
  include '${i}mumps_struc.h'
  external :: ${i}mumps
  type(${i}mumps_struc) :: mumps_par
  end program"
  MUMPS_${i}_links)

if(MUMPS_${i}_links)
  set(MUMPS_${i}_FOUND true PARENT_SCOPE)
  set(MUMPS_links true)
endif()

endforeach()

set(MUMPS_links ${MUMPS_links} PARENT_SCOPE)

endfunction(mumps_check)


function(mumps_libs)

# NOTE: NO_DEFAULT_PATH disables CMP0074 MUMPS_ROOT and PATH_SUFFIXES, so we manually specify:
# HINTS ${MUMPS_ROOT} ENV MUMPS_ROOT
# PATH_SUFFIXES ...
# to allow MKL using user-built MUMPS with `cmake -DMUMPS_ROOT=~/lib_intel/mumps`

if(DEFINED ENV{MKLROOT})
  find_path(MUMPS_INCLUDE_DIR
    NAMES mumps_compat.h
    NO_DEFAULT_PATH
    HINTS ${MUMPS_ROOT} ENV MUMPS_ROOT
    PATH_SUFFIXES include
    DOC "MUMPS common header")
else()
  find_path(MUMPS_INCLUDE_DIR
    NAMES mumps_compat.h
    PATH_SUFFIXES MUMPS openmpi-x86_64 mpich-x86_64
    DOC "MUMPS common header")
endif()
if(NOT MUMPS_INCLUDE_DIR)
  return()
endif()

# --- Mumps Common ---
if(DEFINED ENV{MKLROOT})
  find_library(MUMPS_COMMON
    NAMES mumps_common
    NO_DEFAULT_PATH
    HINTS ${MUMPS_ROOT} ENV MUMPS_ROOT
    PATH_SUFFIXES lib
    DOC "MUMPS MPI common libraries")
elseif(mpiseq IN_LIST MUMPS_FIND_COMPONENTS)
  find_library(MUMPS_COMMON
    NAMES mumps_common mumps_common_seq
    NAMES_PER_DIR
    DOC "MUMPS no-MPI common libraries")
else()
  find_library(MUMPS_COMMON
    NAMES mumps_common mumps_common_mpi mumpso_common mumps_common_shm
    NAMES_PER_DIR
    PATH_SUFFIXES openmpi/lib mpich/lib
    DOC "MUMPS common libraries")
endif()

if(NOT MUMPS_COMMON)
  return()
endif()

# --- Pord ---

if(DEFINED ENV{MKLROOT})
  find_library(PORD
    NAMES pord
    NO_DEFAULT_PATH
    HINTS ${MUMPS_ROOT} ENV MUMPS_ROOT
    PATH_SUFFIXES lib
    DOC "simplest MUMPS ordering library")
else()
  find_library(PORD
    NAMES pord mumps_pord
    NAMES_PER_DIR
    PATH_SUFFIXES openmpi/lib mpich/lib
    DOC "simplest MUMPS ordering library")
endif()
if(NOT PORD)
  return()
endif()

if(mpiseq IN_LIST MUMPS_FIND_COMPONENTS)
  if(DEFINED ENV{MKLROOT})
    find_library(MUMPS_mpiseq_LIB
      NAMES mpiseq
      NO_DEFAULT_PATH
      HINTS ${MUMPS_ROOT} ENV MUMPS_ROOT
      PATH_SUFFIXES lib
      DOC "No-MPI stub library")
  else()
    find_library(MUMPS_mpiseq_LIB
    NAMES mpiseq mumps_mpi_seq
    NAMES_PER_DIR
    DOC "No-MPI stub library")
  endif()
  if(NOT MUMPS_mpiseq_LIB)
    return()
  endif()

  if(DEFINED ENV{MKLROOT})
    find_path(MUMPS_mpiseq_INC
      NAMES mpif.h
      NO_DEFAULT_PATH
      HINTS ${MUMPS_ROOT} ENV MUMPS_ROOT
      PATH_SUFFIXES include
      DOC "MUMPS mpiseq header")
  else()
    find_path(MUMPS_mpiseq_INC
      NAMES mpif.h
      PATH_SUFFIXES MUMPS mumps/mpi_seq
      DOC "MUMPS mpiseq header")
  endif()
  if(NOT MUMPS_mpiseq_INC)
    return()
  endif()

  set(MUMPS_mpiseq_FOUND true PARENT_SCOPE)
  set(MUMPS_mpiseq_LIB ${MUMPS_mpiseq_LIB} PARENT_SCOPE)
  set(MUMPS_mpiseq_INC ${MUMPS_mpiseq_INC} PARENT_SCOPE)
endif()

set(_ariths s d c z)
foreach(comp ${MUMPS_FIND_COMPONENTS})
  if(NOT "${comp}" IN_LIST _ariths)
    continue()
  endif()

  if(DEFINED ENV{MKLROOT})
    find_library(MUMPS_${comp}_lib
      NAMES ${comp}mumps
      NO_DEFAULT_PATH
      HINTS ${MUMPS_ROOT} ENV MUMPS_ROOT
      PATH_SUFFIXES lib
      DOC "MUMPS precision-specific")
  elseif(mpiseq IN_LIST MUMPS_FIND_COMPONENTS)
    find_library(MUMPS_${comp}_lib
      NAMES ${comp}mumps ${comp}mumps_seq
      NAMES_PER_DIR
      DOC "MUMPS no-MPI precision-specific")
  else()
    find_library(MUMPS_${comp}_lib
      NAMES ${comp}mumps ${comp}mumps_mpi
      NAMES_PER_DIR
      PATH_SUFFIXES openmpi/lib mpich/lib
      DOC "MUMPS precision-specific")
  endif()

  if(NOT MUMPS_${comp}_lib)
    return()
  endif()

  set(MUMPS_${comp}_FOUND true PARENT_SCOPE)
  list(APPEND MUMPS_LIBRARY ${MUMPS_${comp}_lib})
endforeach()

set(MUMPS_LIBRARY ${MUMPS_LIBRARY} ${MUMPS_COMMON} ${PORD} PARENT_SCOPE)
set(MUMPS_INCLUDE_DIR ${MUMPS_INCLUDE_DIR} PARENT_SCOPE)

endfunction(mumps_libs)

# --- main

if(NOT MUMPS_FIND_COMPONENTS)
  set(MUMPS_FIND_COMPONENTS d)
endif()

mumps_libs()

if(MUMPS_LIBRARY AND MUMPS_INCLUDE_DIR)

# -- minimal check that MUMPS is linkable

mumps_check()

endif(MUMPS_LIBRARY AND MUMPS_INCLUDE_DIR)
# --- finalize

set(CMAKE_REQUIRED_FLAGS)
set(CMAKE_REQUIRED_INCLUDES)
set(CMAKE_REQUIRED_LIBRARIES)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MUMPS
  REQUIRED_VARS MUMPS_LIBRARY MUMPS_INCLUDE_DIR MUMPS_links
  HANDLE_COMPONENTS)

if(MUMPS_FOUND)
# need if _FOUND guard to allow project to autobuild; can't overwrite imported target even if bad
set(MUMPS_LIBRARIES ${MUMPS_LIBRARY})
set(MUMPS_INCLUDE_DIRS ${MUMPS_INCLUDE_DIR})
if(mpiseq IN_LIST MUMPS_FIND_COMPONENTS)
  list(APPEND MUMPS_LIBRARIES ${MUMPS_mpiseq_LIB})
  list(APPEND MUMPS_INCLUDE_DIRS ${MUMPS_mpiseq_INC})
endif()

if(NOT TARGET MUMPS::MUMPS)
  add_library(MUMPS::MUMPS INTERFACE IMPORTED)
  set_target_properties(MUMPS::MUMPS PROPERTIES
    INTERFACE_LINK_LIBRARIES "${MUMPS_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${MUMPS_INCLUDE_DIR}")
endif()

if(mpiseq IN_LIST MUMPS_FIND_COMPONENTS)
  if(NOT TARGET MUMPS::MPISEQ)
    add_library(MUMPS::MPISEQ INTERFACE IMPORTED)
    set_target_properties(MUMPS::MPISEQ PROPERTIES
      INTERFACE_LINK_LIBRARIES "${MUMPS_mpiseq_LIB}"
      INTERFACE_INCLUDE_DIRECTORIES "${MUMPS_mpiseq_INC}")
  endif()
endif()

endif(MUMPS_FOUND)

mark_as_advanced(MUMPS_INCLUDE_DIR MUMPS_LIBRARY)