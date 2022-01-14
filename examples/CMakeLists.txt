# list(APPEND Examples blasfeo_wrapper_example casadi_codegen_wrapper_example test)
list(APPEND Examples LinearAlgebraExample OCPExample SparseExample OCPBasicExample)
# list(APPEND Examples LinearAlgebraExample OCPExample)

enable_testing()
add_executable(LinearAlgebraExample LinearAlgebraExample.cpp)
add_executable(OCPExample OCPExample.cpp)
add_executable(SparseExample SparseExample.cpp)
add_executable(OCPBasicExample OCPBasicExample.cpp)
# find_package(MUMPS COMPONENTS MUMPS::MPISEQ MUMPS::DMUMPS)
find_package(MUMPS)
target_link_libraries(SparseExample MUMPS::MPISEQ)
target_link_libraries(SparseExample MUMPS::DMUMPS)
# target_link_libraries(SparseExample profiler)


foreach(example IN LISTS Examples) 
    target_link_libraries(${example} fatrop)
    target_include_directories(${example}
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)
endforeach()
