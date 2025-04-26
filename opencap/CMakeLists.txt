cmake_minimum_required(VERSION 3.12)
project(opencap LANGUAGES CXX C)

# Default build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Options
option(BUILD_PYOPENCAP "Builds Python module" OFF)
option(BUILD_OPENCAP "Builds command line version" ON)
option(ENABLE_TESTS "Enable unit tests" ON)

# Language standards and global flags
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Include directories and source files
include_directories(include)
file(GLOB SOURCES "src/*.cpp")
file(GLOB BIND_SOURCES "bindings/*.cpp")

# Dependencies
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})
find_package(OpenMP REQUIRED)

if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Add ASA subdirectory
add_subdirectory(asa)

# External dependencies
include(FetchContent)

FetchContent_Declare(
    h5pp
    GIT_REPOSITORY https://github.com/DavidAce/h5pp.git
    GIT_TAG 0d1cb5ebb6805b3d00fb323a31e91fb5998c7853
)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.14.0
)

FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11
    GIT_TAG v2.13.6
)

FetchContent_Declare(
    numgrid
    GIT_REPOSITORY https://github.com/dftlibs/numgrid.git
    GIT_TAG fdb7b8a0eb482a41b11bfb73b6eff5b0f2611bcf
)

# Disable unit tests for numgrid
set(ENABLE_UNIT_TESTS OFF CACHE INTERNAL "Disable unit tests")
FetchContent_MakeAvailable(h5pp googletest pybind11 numgrid)

# Code coverage (if enabled)
if(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    message(STATUS "Code coverage activated.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -fprofile-arcs -ftest-coverage")
endif()

# RPATH settings
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Build Python module
if(BUILD_PYOPENCAP)
    message(STATUS "Building PyOpenCAP Python module.")
    pybind11_add_module(pyopencap_cpp ${SOURCES} ${BIND_SOURCES})
    target_link_libraries(pyopencap_cpp
        PRIVATE ${HDF5_LIBRARIES} numgrid-objects OpenMP::OpenMP_CXX
        Eigen3::Eigen h5pp::h5pp asa239)
endif()

# Build CLI and shared library
if(BUILD_OPENCAP)
    message(STATUS "Building OpenCAP executable.")
    find_package(Python REQUIRED COMPONENTS Interpreter Development)

    add_executable(opencap ${SOURCES})
    add_library(opencap-shared ${SOURCES})

    target_link_libraries(opencap
        PRIVATE ${HDF5_LIBRARIES} numgrid-objects pybind11::module
        pybind11::embed Eigen3::Eigen h5pp::h5pp OpenMP::OpenMP_CXX asa239)

    target_link_libraries(opencap-shared
        PRIVATE ${HDF5_LIBRARIES} numgrid-objects pybind11::module
        pybind11::embed Eigen3::Eigen h5pp::h5pp OpenMP::OpenMP_CXX asa239)

    install(TARGETS opencap DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

    if(ENABLE_TESTS)
        enable_testing()
        add_subdirectory(tests)
    endif()
endif()
