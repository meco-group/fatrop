# project(spectool)


add_library(spectool SHARED spec/ocp.cpp spec/ustage.cpp solver_interfaces/fatrop/fatrop_solver.cpp solver_interfaces/solver.cpp function_evaluation/casadi_fe.cpp spec/ustage_quantities.cpp spec/stage.cpp spec/ustage_eval_casadi.cpp function_evaluation/casadi_fe.cpp function_evaluation/casadi_jit.cpp)


target_include_directories(spectool PUBLIC  
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>  
    $<INSTALL_INTERFACE:include># <prefix>/include/mylib
)

# # add option for specifying casadi headers location
# set(CASADI_SOURCE "" CACHE STRING "Path to casadi headers")
# # add option for specifying casadi library location
# set(CASADI_LIB "" CACHE STRING "Path to casadi library")
# if(NOT CASADI_SOURCE)
#     message(FATAL_ERROR "Please specify the path to casadi headers")
# endif()
# if(NOT CASADI_LIB)
#     message(FATAL_ERROR "Please specify the path to casadi library")
# endif()
target_link_libraries(spectool fatrop blasfeo casadi)
target_include_directories(spectool PUBLIC ${CASADI_INCLUDE_DIR})
# target_link_libraries(spectool ${CASADI_LIB})
add_executable(example examples/example.cpp)
target_link_libraries(example spectool) 

set(SPECTOOL_INSTALL_HEADERS
auxiliary/casadi_utilities.hpp
auxiliary/constraint_helper.hpp
spec/ocp.hpp
spec/ustage.hpp
solver_interfaces/fatrop/fatrop_solver.hpp
solver_interfaces/solver.hpp
solver_interfaces/fatrop/fatrop_ocp_impl.hpp
solver_interfaces/fatrop/fatrop_function.hpp
spec/ustage_quantities.hpp
spec/stage.hpp
spec/integrator-rk4.hpp
spec/ustage_eval_casadi.hpp
function_evaluation/casadi_fe.hpp
function_evaluation/casadi_jit.hpp
spectool.hpp
)
# install spectool.hpp
install(FILES spectool.hpp DESTINATION include/fatrop/spectool)

# install(TARGETS spectool fatrop blasfeo EXPORT fatropConfig
#     LIBRARY DESTINATION lib
#     ARCHIVE DESTINATION lib
#     RUNTIME DESTINATION bin
#     INCLUDES DESTINATION include
# )

install(DIRECTORY auxiliary spec solver_interfaces function_evaluation DESTINATION include/fatrop/spectool
    FILES_MATCHING PATTERN "*.hpp")

# install(EXPORT spectoolConfig
#     # NAMESPACE fatrop::
#     DESTINATION cmake
# )
