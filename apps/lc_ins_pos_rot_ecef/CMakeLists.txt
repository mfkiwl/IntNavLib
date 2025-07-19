cmake_minimum_required(VERSION 3.7)

project(lc_ins_pos_rot_ecef)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 REQUIRED)
find_package(glog REQUIRED)

find_package(intnavlib REQUIRED)

add_executable(lc_ins_pos_rot_ecef 
                src/lc_ins_pos_rot_ecef.cpp
                )

target_link_libraries(lc_ins_pos_rot_ecef  
                                    PRIVATE
                                        ${OPENTHREADS_LIBRARY}
                                        Eigen3::Eigen
                                        glog::glog
                                        intnavlib::intnavlib
                                        )