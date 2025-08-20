cmake_minimum_required(VERSION 3.7)

project(nav_sim)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 REQUIRED)
find_package(glog REQUIRED)

find_package(intnavlib REQUIRED)

add_executable(nav_sim 
                    src/nav_sim.cpp
                    )

target_link_libraries(nav_sim  
                                    PRIVATE
                                        ${OPENTHREADS_LIBRARY}
                                        Eigen3::Eigen
                                        glog::glog
                                        intnavlib::intnavlib
                                        )