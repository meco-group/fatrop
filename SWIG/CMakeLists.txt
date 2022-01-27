add_custom_command(OUTPUT BenchmarkAux.cxx COMMAND swig -o BenchmarkAux.cxx -python -c++ ${CMAKE_SOURCE_DIR}/SWIG/BenchmarkAux.i  DEPENDS BenchmarkAux.i)
add_library(_BenchmarkAux SHARED BenchmarkAux.cxx)
target_link_libraries(_BenchmarkAux fatrop dl) 
find_package (Python3 COMPONENTS Interpreter Development NumPy)
# find_package (Python3 COMPONENTS Interpreter Development) 
message("Python_FOUND:${Python3_FOUND}")
message("Python_VERSION:${Python3_VERSION}")
message("Python_Development_FOUND:${Python3_Development_FOUND}")
message("Python_LIBRARIES:${Python3_LIBRARIES}")
target_include_directories(_BenchmarkAux
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${Python3_INCLUDE_DIRS}
        ${Python3_NumPy_INCLUDE_DIRS}
)
set_target_properties(_BenchmarkAux PROPERTIES PREFIX "")