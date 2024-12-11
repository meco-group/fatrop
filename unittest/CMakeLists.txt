find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_executable(FatropOptionsTests FatropOptionsTests.cpp)
target_link_libraries(FatropOptionsTests fatrop)
target_link_libraries(FatropOptionsTests ${GTEST_LIBRARIES})