add_library(spectool SHARED ocp.cpp ustage.cpp fatrop_solver.cpp solver.cpp casadi_fe.cpp fatrop_ocp_impl.cpp ustage_quantities.cpp stage.cpp)
target_link_libraries(spectool fatrop blasfeo casadi)
add_executable(example example.cpp)
target_link_libraries(example spectool)