add_executable(matrix_vector_example matrix_vector_example.cpp)
target_link_libraries(matrix_vector_example PRIVATE blasfeo)

add_executable(ocp_impl_example ocp_impl_example.cpp)
target_link_libraries(ocp_impl_example PRIVATE fatrop)

add_executable(ocp_impl_example_static ocp_impl_example_static.cpp)
target_link_libraries(ocp_impl_example_static PRIVATE fatrop)

# Include the casadi_codegen subdirectory
add_subdirectory(casadi_codegen)
