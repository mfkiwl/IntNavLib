cmake_minimum_required(VERSION 3.7)

project(tc_ins_gnss_ecef)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 REQUIRED)
find_package(glog REQUIRED)

find_package(intnavlib REQUIRED)

add_executable(tc_ins_gnss_ecef 
                    src/tc_ins_gnss_ecef.cpp
                    )

target_link_libraries(tc_ins_gnss_ecef  
                                    PRIVATE
                                        ${OPENTHREADS_LIBRARY}
                                        Eigen3::Eigen
                                        glog::glog
                                        intnavlib::intnavlib
                                        )