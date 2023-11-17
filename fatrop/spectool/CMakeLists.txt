add_library(spectool SHARED spec/ocp.cpp spec/ustage.cpp solver_interfaces/fatrop/fatrop_solver.cpp solver_interfaces/solver.cpp function_evaluation/casadi_fe.cpp solver_interfaces/fatrop/fatrop_ocp_impl.cpp spec/ustage_quantities.cpp spec/stage.cpp solver_interfaces/fatrop/fatrop_ustage_eval_casadi.cpp function_evaluation/casadi_fe.cpp function_evaluation/casadi_jit.cpp)


target_include_directories(spectool PUBLIC  
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>  
    $<INSTALL_INTERFACE:include> # <prefix>/include/mylib
)

# add option for specifying casadi headers location
set(CASADI_SOURCE "" CACHE STRING "Path to casadi headers")
# add option for specifying casadi library location
set(CASADI_LIB "" CACHE STRING "Path to casadi library")
if(NOT CASADI_SOURCE)
    message(FATAL_ERROR "Please specify the path to casadi headers")
endif()
if(NOT CASADI_LIB)
    message(FATAL_ERROR "Please specify the path to casadi library")
endif()
target_link_libraries(spectool fatrop blasfeo)
target_include_directories(spectool PUBLIC ${CASADI_SOURCE})
target_link_libraries(spectool ${CASADI_LIB})
add_executable(example examples/example.cpp)
target_link_libraries(example spectool)