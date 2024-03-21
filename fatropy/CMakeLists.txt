
pybind11_add_module(_fatropy src/fatropy.cpp)
target_link_libraries(_fatropy PRIVATE fatrop)
if(${WITH_SPECTOOL})
    target_link_libraries(_fatropy PRIVATE spectool casadi)
endif()
target_include_directories(_fatropy PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../fatrop)
target_include_directories(_fatropy PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../fatrop/spectool)

# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
target_compile_definitions(_fatropy
                           PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
# if spectool 
if(${WITH_SPECTOOL})
add_library(python_func_import SHARED src/fatropy.cpp src/python_func_import.cpp)
target_link_libraries(python_func_import PRIVATE fatrop pybind11::module pybind11::embed)
endif()