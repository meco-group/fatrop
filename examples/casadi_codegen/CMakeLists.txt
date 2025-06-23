# Attempt to run the generate_casadi.py script
execute_process(
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/generate_casadi.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    RESULT_VARIABLE GENERATE_CASADI_RESULT
    OUTPUT_VARIABLE GENERATE_CASADI_OUTPUT
    ERROR_VARIABLE GENERATE_CASADI_ERROR
)

# Check if the script ran successfully
if(NOT GENERATE_CASADI_RESULT EQUAL 0)
    message(STATUS "Failed to generate CasADi files. Error: ${GENERATE_CASADI_ERROR}")
    message(STATUS "Skipping casadi code generation example.")
else()
    message(STATUS "Successfully generated CasADi files.")
    
    # Compile the casadi generated c code to a shared library
    add_library(casadi_generated SHARED casadi_generated.c)
    target_link_libraries(casadi_generated PRIVATE fatrop)

    add_executable(ocp_impl_example_codegen ocp_impl_example_codegen.cpp)
    target_link_libraries(ocp_impl_example_codegen PRIVATE casadi_generated fatrop)
    target_include_directories(ocp_impl_example_codegen PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
endif()
