add_executable(graph_example graph_example.cpp)
target_link_libraries(graph_example PRIVATE fatrop)

add_executable(graph_slam_example graph_slam_example.cpp)
target_link_libraries(graph_slam_example PRIVATE fatrop)
target_include_directories(graph_slam_example PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)
