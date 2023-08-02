find_package(pybind11 REQUIRED)

pybind11_add_module(fatropy fatropy.cpp)
target_link_libraries(fatropy PRIVATE fatrop)
target_include_directories(fatropy PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../fatrop)

# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
target_compile_definitions(fatropy
                           PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})