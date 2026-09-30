option(BENCHMARKS_WITH_IPOPT "Build the Ipopt comparison benchmark" OFF)

add_executable(fatrop_ocp_derivative_generator codegen/generate_derivatives.cpp)
target_link_libraries(fatrop_ocp_derivative_generator PRIVATE fatrop)

set(generated_derivatives "${CMAKE_CURRENT_BINARY_DIR}/compiled_derivatives.hpp")
add_custom_command(
    OUTPUT "${generated_derivatives}"
    COMMAND fatrop_ocp_derivative_generator "${generated_derivatives}"
    DEPENDS fatrop_ocp_derivative_generator
    VERBATIM
)

add_executable(fatrop_ocp_bfgs_benchmark
    ocp_bfgs_scenarios.cpp
    "${generated_derivatives}"
)
target_include_directories(fatrop_ocp_bfgs_benchmark PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")
target_link_libraries(fatrop_ocp_bfgs_benchmark PRIVATE fatrop)

add_executable(fatrop_ocp_derivative_check
    codegen/check_compiled_derivatives.cpp
    "${generated_derivatives}"
)
target_include_directories(fatrop_ocp_derivative_check PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")
target_link_libraries(fatrop_ocp_derivative_check PRIVATE fatrop)

if(BENCHMARKS_WITH_IPOPT)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(IPOPT REQUIRED IMPORTED_TARGET ipopt)
    add_executable(fatrop_ocp_vs_ipopt_benchmark ocp_vs_ipopt_benchmark.cpp)
    target_link_libraries(fatrop_ocp_vs_ipopt_benchmark PRIVATE fatrop PkgConfig::IPOPT)
endif()
