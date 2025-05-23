cmake_minimum_required(VERSION 3.7)

project(ins_ecef)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(Eigen3 REQUIRED)
find_package(glog REQUIRED)

find_package(intnavlib REQUIRED)

add_executable(ins_ecef 
                    src/ins_ecef.cpp
                    )

target_include_directories(ins_ecef 
                                    PRIVATE 
                                        ${INTNAVLIB_INCLUDE_DIRS})

target_link_libraries(ins_ecef  
                                    PRIVATE
                                        ${OPENTHREADS_LIBRARY}
                                        Eigen3::Eigen
                                        glog::glog
                                        ${INTNAVLIB_LIBRARIES}
                                        )