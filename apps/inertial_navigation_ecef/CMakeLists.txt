cmake_minimum_required(VERSION 3.7)

project(inertial_navigation_ecef)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(Eigen3 REQUIRED)
find_package(glog REQUIRED)

find_package(intnavlib REQUIRED)

add_executable(inertial_navigation_ecef 
                    src/inertial_navigation_ecef.cpp
                    )

target_include_directories(inertial_navigation_ecef 
                                    PRIVATE 
                                        ${INTNAVLIB_INCLUDE_DIRS})

target_link_libraries(inertial_navigation_ecef  
                                    PRIVATE
                                        ${OPENTHREADS_LIBRARY}
                                        Eigen3::Eigen
                                        glog::glog
                                        ${INTNAVLIB_LIBRARIES}
                                        )