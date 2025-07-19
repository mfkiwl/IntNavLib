cmake_minimum_required(VERSION 3.7)

project(intnavlib VERSION 1.0 DESCRIPTION "integrated navigation c++ library")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Export compile commands for clang tidy
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Eigen3 REQUIRED)

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
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include/intnavlib>
)

# Set properties
set_target_properties(intnavlib PROPERTIES
                        VERSION ${PROJECT_VERSION}
                        SOVERSION ${PROJECT_VERSION_MAJOR}
                        POSITION_INDEPENDENT_CODE ON
                        LINKER_LANGUAGE CXX 
                        PUBLIC_HEADER "${PUBLIC_HEADERS}")


# Link dependencies. PRIVATE so only visible from lib
target_link_libraries( intnavlib
        PRIVATE
            Eigen3::Eigen
        )

# Install instructions
include(GNUInstallDirs)
install(TARGETS intnavlib
    EXPORT intnavlibTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/intnavlib)

# Export targets
install(EXPORT intnavlibTargets
    FILE intnavlibTargets.cmake
    NAMESPACE intnavlib::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/intnavlib
)

# Config files
include(CMakePackageConfigHelpers)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/intnavlibConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/intnavlibConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/intnavlibConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/intnavlib
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/intnavlibConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/intnavlibConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/intnavlib
)
