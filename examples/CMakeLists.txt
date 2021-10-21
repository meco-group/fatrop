# list(APPEND Examples blasfeo_wrapper_example casadi_codegen_wrapper_example test)
list(APPEND Examples LinearAlgebraExample)

add_executable(LinearAlgebraExample LinearAlgebraExample.cpp)

# link_directories(fatrop_core, blasfeo_wrapper)

foreach(example IN LISTS Examples) 
    target_link_libraries(${example} fatrop)
    target_include_directories(${example}
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)
endforeach()