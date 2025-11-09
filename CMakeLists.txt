cmake_minimum_required(VERSION 3.7)

project(intnavlib VERSION 1.0 DESCRIPTION "integrated navigation c++ library")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Export compile commands for clang tidy
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(USE_FLOAT "Build with float precision instead of double" OFF)

if(USE_FLOAT)
    message(STATUS "Build with float precision")
    add_compile_definitions(USE_FLOAT)
endif()

option(CODE_COVERAGE "Enable coverage reporting" OFF)

function(enable_coverage target)
    if(CHECK_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${target} PRIVATE -O0 -g --coverage)
        target_link_options(${target} PRIVATE --coverage)
    endif()
endfunction()

find_package(Eigen3 REQUIRED)

# Get all header files in the include directory
file(GLOB PUBLIC_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h")
# Get all source files in the src directory
file(GLOB PRIVATE_SOURCE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

# Create static library with its sources
# add_library(intnavlib STATIC ${PRIVATE_SOURCE_FILES}) 

add_library(intnavlib INTERFACE) # header-only lib

# Add headers. We choose to keep them all public
# target_include_directories(intnavlib
#     PUBLIC
#         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
#         $<INSTALL_INTERFACE:include/intnavlib>
# )
target_include_directories(intnavlib
    INTERFACE
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
# target_link_libraries( intnavlib
#         PRIVATE
#             Eigen3::Eigen
#         )
target_link_libraries( intnavlib
        INTERFACE
            Eigen3::Eigen
        )


# ================= Tests =================

# https://github.com/google/googletest/tree/main/googletest#incorporating-into-an-existing-cmake-project

find_package(glog REQUIRED)

option(BUILD_TESTING "Enable testing" ON)
if(BUILD_TESTING)
    include(CTest)
    enable_testing()

    # Fetch GoogleTest
    include(FetchContent)
    FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/5376968f6948923e2411081fd9372e71a59d8e77.zip
    )

    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)

    # Get all test sources
    file(GLOB TEST_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp")

    # Creeate tests
    foreach(test_file ${TEST_SOURCES})
        get_filename_component(test_name ${test_file} NAME_WE)
        add_executable(${test_name} ${test_file})
        target_link_libraries(${test_name}
            PRIVATE
                gtest
                gtest_main
                ${OPENTHREADS_LIBRARY}
                Eigen3::Eigen
                glog::glog
                intnavlib
        )
        add_test(NAME ${test_name} COMMAND ${test_name})
        enable_coverage(${test_name})
    endforeach()
endif()

# ================= Install =================

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