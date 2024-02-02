
pybind11_add_module(fatropy fatropy.cpp)
target_link_libraries(fatropy PRIVATE fatrop)
if(${WITH_SPECTOOL})
    target_link_libraries(fatropy PRIVATE spectool casadi)
endif()
target_include_directories(fatropy PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../fatrop)
target_include_directories(fatropy PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../fatrop/spectool)

# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
target_compile_definitions(fatropy
                           PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
add_library(python_func_import SHARED fatropy.cpp src/python_func_import.cpp)
target_link_libraries(python_func_import PRIVATE fatrop pybind11::module pybind11::embed)