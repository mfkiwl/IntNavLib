cmake_minimum_required(VERSION 3.7)

project(intnavlib VERSION 1.0 DESCRIPTION "integrated navigation c++ library")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 REQUIRED)
find_package(glog REQUIRED)

# Create static library
add_library(intnavlib STATIC) 

# Get all header files in the include directory
file(GLOB PUBLIC_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h")

# Get all source files in the src directory
file(GLOB PRIVATE_SOURCE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

# Add private source files
target_sources(intnavlib
    PRIVATE
    ${PRIVATE_SOURCE_FILES}
)

# Add headers. We choose to keep them all public
target_include_directories(intnavlib
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Set properties
set_target_properties(intnavlib PROPERTIES
                        VERSION ${PROJECT_VERSION}
                        SOVERSION ${PROJECT_VERSION_MAJOR}
                        LINKER_LANGUAGE CXX 
                        PUBLIC_HEADER "${PUBLIC_HEADERS}")


# Link dependencies. PRIVATE so only visible from lib
target_link_libraries( intnavlib
        PRIVATE
            Eigen3::Eigen
            glog::glog
        )

# Install instructions
include(GNUInstallDirs)
install(TARGETS intnavlib
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/intnavlib
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/intnavlib)