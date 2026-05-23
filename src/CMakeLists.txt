# Define source files for the fatrop library
set(FATROP_SOURCES
    linear_algebra/lu_factorization.cpp
    linear_algebra/vector.cpp
    ocp/aug_system_solver.cpp
    ocp/dims.cpp
    ocp/jacobian.cpp
    ocp/hessian.cpp
    ocp/ip_alg_builder.cpp
    ocp/ip_algorithm.cpp
    ocp/ip_convergence_check.cpp
    ocp/ip_convergence_check_resto.cpp
    ocp/ip_data.cpp
    ocp/ip_eq_mult_initializer.cpp
    ocp/ip_initializer.cpp
    ocp/ip_initializer_resto.cpp
    ocp/ip_iterate.cpp
    ocp/ip_iteration_output.cpp
    ocp/ip_iteration_output_resto.cpp
    ocp/ip_linesearch.cpp
    ocp/ip_nlp_orig.cpp
    ocp/ip_nlp_resto.cpp
    ocp/ip_resto_phase_min_cl1.cpp
    ocp/ip_search_dir.cpp
    ocp/nlp_ocp.cpp
    ocp/pd_system_orig.cpp
    ocp/pd_system_resto.cpp
    ocp/pd_solver_orig.cpp
    ocp/pd_solver_resto.cpp
    ocp/problem_info.cpp
    ocp/mehrotra_qp.cpp
    ip_algorithm/ip_filter.cpp
    ip_algorithm/ip_mu_update.cpp
    ip_algorithm/ip_utils.cpp
    dense/dims.cpp
    dense/problem_info.cpp
    dense/jacobian.cpp
    dense/hessian.cpp
    dense/aug_system_solver.cpp
    dense/pd_system_orig.cpp
    dense/pd_solver_orig.cpp
    dense/pd_system_resto.cpp
    dense/pd_solver_resto.cpp
    dense/nlp_dense.cpp
    dense/ip_alg_builder.cpp
    dense/ip_algorithm.cpp
    dense/ip_convergence_check.cpp
    dense/ip_convergence_check_resto.cpp
    dense/ip_data.cpp
    dense/ip_eq_mult_initializer.cpp
    dense/ip_initializer.cpp
    dense/ip_initializer_resto.cpp
    dense/ip_iterate.cpp
    dense/ip_iteration_output.cpp
    dense/ip_iteration_output_resto.cpp
    dense/ip_linesearch.cpp
    dense/ip_mu_update.cpp
    dense/ip_nlp_orig.cpp
    dense/ip_nlp_resto.cpp
    dense/ip_resto_phase_min_cl1.cpp
    dense/ip_search_dir.cpp
    dense/mehrotra_qp.cpp
    graph/block_sparsity.cpp
    graph/block_pd_matrix.cpp
    graph/linear_system.cpp
    graph/elimination_order.cpp
    graph/symbolic_factorization.cpp
    graph/linear_solver.cpp
    graph/dims.cpp
    graph/problem_info.cpp
    graph/hessian.cpp
    graph/jacobian.cpp
    graph/nlp_graph.cpp
    graph/aug_system_solver.cpp
    graph/pd_system_orig.cpp
    graph/pd_solver_orig.cpp
    graph/pd_system_resto.cpp
    graph/pd_solver_resto.cpp
    graph/ip_iterate.cpp
    graph/ip_data.cpp
    graph/ip_nlp_orig.cpp
    graph/ip_nlp_resto.cpp
    graph/ip_search_dir.cpp
    graph/ip_linesearch.cpp
    graph/ip_mu_update.cpp
    graph/ip_convergence_check.cpp
    graph/ip_convergence_check_resto.cpp
    graph/ip_iteration_output.cpp
    graph/ip_iteration_output_resto.cpp
    graph/ip_eq_mult_initializer.cpp
    graph/ip_initializer.cpp
    graph/ip_initializer_resto.cpp
    graph/ip_resto_phase_min_cl1.cpp
    graph/ip_algorithm.cpp
    graph/ip_alg_builder.cpp
)

# Create the fatrop shared library
add_library(fatrop SHARED ${FATROP_SOURCES})

# Set properties for the fatrop target
set_target_properties(fatrop PROPERTIES
    VERSION ${PROJECT_VERSION}
)

# Link fatrop with blasfeo
target_link_libraries(fatrop PUBLIC blasfeo)

# Specify include directories for the target
target_include_directories(fatrop
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)


# Define install targets
set(INSTALL_TARGETS fatrop)

# Export the target
install(TARGETS ${INSTALL_TARGETS}
    EXPORT fatropConfig
    LIBRARY DESTINATION ${LIBRARY_INSTALL_DIR}
    ARCHIVE DESTINATION ${ARCHIVE_INSTALL_DIR}
    RUNTIME DESTINATION ${RUNTIME_INSTALL_DIR}
    INCLUDES DESTINATION ${INCLUDE_INSTALL_DIR}
)

# Install header files (only .hpp)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/fatrop
        DESTINATION ${INCLUDE_INSTALL_DIR}
        FILES_MATCHING PATTERN "*.hpp")

# Install the export
install(EXPORT fatropConfig
        NAMESPACE fatrop::
        DESTINATION ${CMAKE_INSTALL_DIR})
