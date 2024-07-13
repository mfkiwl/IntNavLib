cmake_minimum_required(VERSION 3.7)

project(lc_ins_gnss)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(Eigen3 REQUIRED)
find_package(glog REQUIRED)

find_package(intnavlib REQUIRED)

add_executable(lc_ins_gnss 
                    src/lc_ins_gnss.cpp
                    )

target_include_directories(lc_ins_gnss 
                                    PRIVATE 
                                        ${INTNAVLIB_INCLUDE_DIRS})

target_link_libraries(lc_ins_gnss  
                                    PRIVATE
                                        ${OPENTHREADS_LIBRARY}
                                        Eigen3::Eigen
                                        glog::glog
                                        ${INTNAVLIB_LIBRARIES}
                                        )