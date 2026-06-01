add_executable(dense_example dense_example.cpp)
target_link_libraries(dense_example PRIVATE fatrop)

add_executable(quaternion_lie_example quaternion_lie_example.cpp)
target_link_libraries(quaternion_lie_example PRIVATE fatrop)
target_include_directories(quaternion_lie_example PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)

add_executable(pnp_lie_example pnp_lie_example.cpp)
target_link_libraries(pnp_lie_example PRIVATE fatrop)
target_include_directories(pnp_lie_example PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)

add_executable(robot_ik_example robot_ik_example.cpp)
target_link_libraries(robot_ik_example PRIVATE fatrop)
target_include_directories(robot_ik_example PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)
