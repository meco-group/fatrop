# list(APPEND Examples blasfeo_wrapper_example casadi_codegen_wrapper_example test)
list(APPEND Tests LinearAlgebraTests)

enable_testing()
add_executable(LinearAlgebraTests LinearAlgebraTests.cpp)
target_link_libraries(
  LinearAlgebraTests
  gtest_main
)

foreach(Test IN LISTS Tests) 
    target_link_libraries(${Test} fatrop)
    target_include_directories(${Test}
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)
endforeach()

include(GoogleTest)
gtest_discover_tests(LinearAlgebraTests)
