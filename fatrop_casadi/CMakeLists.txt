add_library(fatrop_casadi SHARED specification.cpp)
target_link_libraries(fatrop_casadi fatrop blasfeo casadi)
add_executable(example example.cpp)
target_link_libraries(example fatrop_casadi)