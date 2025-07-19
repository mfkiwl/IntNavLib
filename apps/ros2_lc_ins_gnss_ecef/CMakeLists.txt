cmake_minimum_required(VERSION 3.7)

project(ros2_lc_ins_gnss_ecef)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall) # Wextra -Wpedantic
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  add_compile_options(-O0) # Do not optimize
else()
  add_compile_options(-O3) # Optimize
endif()

# ROS2 dependencies
find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(std_msgs REQUIRED)
find_package(sensor_msgs REQUIRED)
find_package(message_filters REQUIRED)
find_package(visualization_msgs REQUIRED)
find_package(nav_msgs REQUIRED)
find_package(tf2 REQUIRED)
find_package(tf2_ros REQUIRED)
find_package(geometry_msgs REQUIRED)

find_package(Eigen3 REQUIRED)
find_package(glog REQUIRED)

find_package(intnavlib REQUIRED)

add_executable(ros2_lc_ins_gnss_ecef 
                    src/ros2_lc_ins_gnss_ecef.cpp
                    )

target_link_libraries(ros2_lc_ins_gnss_ecef  
                                        ${OPENTHREADS_LIBRARY}
                                        Eigen3::Eigen
                                        glog::glog
                                        intnavlib::intnavlib
                                        )

ament_target_dependencies(ros2_lc_ins_gnss_ecef 
                            rclcpp 
                            std_msgs 
                            sensor_msgs
                            nav_msgs
                            geometry_msgs
                            message_filters
                            visualization_msgs
                            tf2
                            tf2_ros)

############### install ##############

# This puts include files, libs and launch files in the install dir in the ros ws

# Install launch files.
install(DIRECTORY
  launch
  DESTINATION share/${PROJECT_NAME}/
)

# Install libraries
install(TARGETS ros2_lc_ins_gnss_ecef
  DESTINATION lib/${PROJECT_NAME}/
)

# # Install include directories
# install(DIRECTORY include/${PROJECT_NAME}/
#   DESTINATION include/${PROJECT_NAME}
# )

# Install python programs
# install(PROGRAMS
#   ${PROJECT_NAME}/whole_body_robot_bringup.py
#   DESTINATION lib/${PROJECT_NAME}
# )

ament_package()