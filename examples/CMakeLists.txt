# list(APPEND Examples blasfeo_wrapper_example casadi_codegen_wrapper_example test)
list(APPEND Examples LinearAlgebraExample OCPExample)

enable_testing()
add_executable(LinearAlgebraExample LinearAlgebraExample.cpp)
add_executable(OCPExample OCPExample.cpp)

foreach(example IN LISTS Examples) 
    target_link_libraries(${example} fatrop)
    target_include_directories(${example}
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)
endforeach()
