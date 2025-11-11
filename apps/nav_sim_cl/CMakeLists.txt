cmake_minimum_required(VERSION 3.7)

project(nav_sim_cl)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 REQUIRED)
find_package(glog REQUIRED)

find_package(intnavlib REQUIRED)

add_executable(nav_sim_cl 
                    src/nav_sim_cl.cpp
                    )

target_link_libraries(nav_sim_cl  
                        PRIVATE
                            ${OPENTHREADS_LIBRARY}
                            Eigen3::Eigen
                            glog::glog
                            intnavlib::intnavlib
                            )

# manually include and link JSBSim
target_include_directories(nav_sim_cl PRIVATE /usr/include/JSBSim)
target_link_libraries(nav_sim_cl PRIVATE /usr/lib/libJSBSim.a)