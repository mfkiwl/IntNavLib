cmake_minimum_required(VERSION 3.7)

project(intnavlib VERSION 1.0 DESCRIPTION "integrated navigation c++ library")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 REQUIRED)
find_package (glog REQUIRED)

# Create static library
add_library(intnavlib STATIC) 

# Add private source files
# private sources are only visible inside lib
target_sources(intnavlib
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Add headers. We choose to keep them all public
# public means visible from outside
target_include_directories(intnavlib
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Set properties
set_target_properties(intnavlib PROPERTIES
                        VERSION ${PROJECT_VERSION}
                        SOVERSION ${PROJECT_VERSION_MAJOR}
                        LINKER_LANGUAGE CXX 
                        PUBLIC_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/include/intnavlib.h)


# Link dependencies. PRIVATE so only visible from lib
target_link_libraries( intnavlib
        PRIVATE
            Eigen3
            glog
        )

# Install instructions
include(GNUInstallDirs)
install(TARGETS intnavlib
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Optionally to install files such as config files, readme, ...
# install(FILES ${CMAKE_BINARY_DIR}/mylib.pc
#     DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)