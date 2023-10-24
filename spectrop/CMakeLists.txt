add_library(spectrop SHARED ocp.cpp stage.cpp fatrop_solver.cpp solver.cpp casadi_fe.cpp fatrop_ocp_impl.cpp stage_quantities.cpp)
target_link_libraries(spectrop fatrop blasfeo casadi)
add_executable(example example.cpp)
target_link_libraries(example spectrop)