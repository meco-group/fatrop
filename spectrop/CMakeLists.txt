add_library(spectrop SHARED ocp.cpp stage.cpp solver_fatrop.cpp solver.cpp)
target_link_libraries(spectrop fatrop blasfeo casadi)
add_executable(example example.cpp)
target_link_libraries(example spectrop)