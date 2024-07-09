cmake_minimum_required(VERSION 3.7)
project(intnavlib)

# set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 REQUIRED)

find_package (glog REQUIRED)

# include the include folder
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Create library


# ======== target 1: ins_ecef_demo ========

add_executable(ins_ecef_demo 
                    src/ins_ecef_demo.cpp
                    src/helpers.cpp
                    src/lc_ins_gnss
                    )

target_link_libraries(ins_ecef_demo PRIVATE  ${OPENTHREADS_LIBRARY}
                                    ${OpenCV_LIBS}
                                    Eigen3::Eigen
                                    glog::glog
                                    )