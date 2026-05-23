add_executable(matrix_vector_example matrix_vector_example.cpp)
target_link_libraries(matrix_vector_example PRIVATE blasfeo)

add_subdirectory(ocp)
add_subdirectory(dense)
add_subdirectory(graph)
add_subdirectory(casadi)
