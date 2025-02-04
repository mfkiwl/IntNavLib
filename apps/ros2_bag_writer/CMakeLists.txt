cmake_minimum_required(VERSION 3.8)
project(ros2_bag_writer)

if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
endif()

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(rosbag2_cpp REQUIRED)
find_package(geometry_msgs REQUIRED)
find_package(tf2_msgs REQUIRED)
find_package(sensor_msgs REQUIRED)

find_package(Eigen3 REQUIRED)
find_package(glog REQUIRED)
find_package(intnavlib REQUIRED)

add_executable(ros2_bag_writer src/ros2_bag_writer.cpp)

target_include_directories(ros2_bag_writer 
                                    PRIVATE 
                                        ${INTNAVLIB_INCLUDE_DIRS})

target_link_libraries(ros2_bag_writer  
                      Eigen3::Eigen
                      glog::glog
                      ${INTNAVLIB_LIBRARIES})

ament_target_dependencies(ros2_bag_writer
  rclcpp
  rosbag2_cpp
  geometry_msgs
  tf2_msgs
  sensor_msgs
)

install(TARGETS ros2_bag_writer
  DESTINATION lib/${PROJECT_NAME})


ament_package()