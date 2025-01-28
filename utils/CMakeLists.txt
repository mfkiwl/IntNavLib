cmake_minimum_required(VERSION 3.7)

project(utils)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(Eigen3 REQUIRED)
find_package(glog REQUIRED)

find_package(intnavlib REQUIRED)

add_executable(compute_errors 
                    compute_errors.cpp
                    )

target_include_directories(compute_errors 
                                    PRIVATE 
                                        ${INTNAVLIB_INCLUDE_DIRS})

target_link_libraries(compute_errors  
                                    PRIVATE
                                        Eigen3::Eigen
                                        ${INTNAVLIB_LIBRARIES}
                                        )