add_library(spectrop SHARED ocp.cpp ustage.cpp fatrop_solver.cpp solver.cpp casadi_fe.cpp fatrop_ocp_impl.cpp ustage_quantities.cpp stage.cpp)
target_link_libraries(spectrop fatrop blasfeo casadi)
add_executable(example example.cpp)
target_link_libraries(example spectrop)