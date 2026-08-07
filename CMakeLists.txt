cmake_minimum_required(VERSION 3.14)
project(fatrop VERSION 1.1.6 LANGUAGES CXX C)

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

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

# Set compilation flags for Debug
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")

# Option to build tests (default is ON)
option(BUILD_TESTS "Build the tests" ON)
option(BUILD_EXAMPLES "Build the examples" ON)
option(BUILD_WITH_C_INTERFACE "Build with C interface support" ON)
option(FATROP_ENABLE_ASAN "Build with AddressSanitizer/UBSan and force-enable assertions" OFF)

if(FATROP_ENABLE_ASAN)
  add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer -g -UNDEBUG)
  add_link_options(-fsanitize=address,undefined)
endif()

# Define installation directories
set(LIBRARY_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}")
set(ARCHIVE_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}")
set(RUNTIME_INSTALL_DIR "${CMAKE_INSTALL_BINDIR}")
set(INCLUDE_INSTALL_DIR "${CMAKE_INSTALL_INCLUDEDIR}")
set(CMAKE_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/fatropConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/fatropConfig.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_DIR}"
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/fatropConfigVersion.cmake"
    VERSION "${PROJECT_VERSION}"
    COMPATIBILITY SameMajorVersion
)

install(
    FILES
        "${CMAKE_CURRENT_BINARY_DIR}/fatropConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/fatropConfigVersion.cmake"
    DESTINATION "${CMAKE_INSTALL_DIR}"
)

# Select the linear-algebra backend
set(FATROP_LA_BACKEND "blasfeo" CACHE STRING "Linear algebra backend (blasfeo or eigen)")
set_property(CACHE FATROP_LA_BACKEND PROPERTY STRINGS "blasfeo" "eigen")

# Add the external directory
add_subdirectory(external)

if(FATROP_LA_BACKEND STREQUAL "eigen")
    message(STATUS "Building fatrop with the Eigen linear-algebra backend")
    find_package(Eigen3 3.3 QUIET NO_MODULE)
    if(NOT TARGET Eigen3::Eigen)
        message(STATUS "System Eigen3 not found. Fetching Eigen 3.4.0.")
        include(FetchContent)
        set(EIGEN_BUILD_DOC OFF)
        set(EIGEN_BUILD_TESTING OFF)
        set(EIGEN_BUILD_PKGCONFIG OFF)
        FetchContent_Declare(
            eigen
            GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
            GIT_TAG 3.4.0
        )
        FetchContent_MakeAvailable(eigen)
    endif()
    if(BUILD_WITH_C_INTERFACE)
        message(STATUS "The C interface exposes blasfeo types in its ABI and is only available "
                       "with the blasfeo backend - disabling it")
        set(BUILD_WITH_C_INTERFACE OFF)
    endif()
else()
    # Add BLASFEO
    if(WITH_BUILD_BLASFEO)
        # BLASFEO is already available from FetchContent
        add_library(blasfeo::blasfeo ALIAS blasfeo)
    else()
        find_package(blasfeo REQUIRED)
    endif()
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
