cmake_minimum_required(VERSION 3.12)
project(opencap
        LANGUAGES CXX C 
        )
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
option(BUILD_PYOPENCAP    "Builds Python module"          OFF)
option(BUILD_OPENCAP      "Builds command line version"   ON)
include_directories(include)
file(GLOB SOURCES "src/*.cpp")
file(GLOB BIND_SOURCES "bindings/*.cpp")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})
include(FetchContent)


find_package(OpenMP REQUIRED)

FetchContent_Declare(
    h5pp
    GIT_REPOSITORY https://github.com/DavidAce/h5pp.git
    GIT_TAG v1.7.3
)
FetchContent_GetProperties(h5pp)
if(NOT h5pp_POPULATED)
    FetchContent_Populate(h5pp)
    add_subdirectory(${h5pp_SOURCE_DIR} ${h5pp_BINARY_DIR})
endif()

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        703bd9caab50b139428cea1aaff9974ebee5742e # release-1.10.0
)
if(NOT googletest_POPULATED)
  FetchContent_Populate(googletest)
  add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})
endif()


FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11
    GIT_TAG v2.6.2
)
FetchContent_GetProperties(pybind11)
if(NOT pybind11_POPULATED)
    FetchContent_Populate(pybind11)
    add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
endif()

FetchContent_Declare(
  numgrid
  GIT_REPOSITORY https://github.com/dftlibs/numgrid.git
  GIT_TAG fdb7b8a0eb482a41b11bfb73b6eff5b0f2611bcf
)
FetchContent_GetProperties(numgrid)
if(NOT numgrid_POPULATED)
  FetchContent_Populate(numgrid)
  set(ENABLE_UNIT_TESTS OFF CACHE INTERNAL "Enable unit tests")
  add_subdirectory(${numgrid_SOURCE_DIR} ${numgrid_BINARY_DIR})
endif()

# Code Coverage Configuration
if(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    MESSAGE("Code coverage activated.")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -fprofile-arcs -ftest-coverage")
	SET(CMAKE_CXX_OUTPUT_EXTENSION_REPLACE 1)
endif(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

if(BUILD_PYOPENCAP)
message("Building PyOpenCAP python module.")
pybind11_add_module(pyopencap_cpp ${SOURCES} ${BIND_SOURCES})
target_link_libraries(pyopencap_cpp PRIVATE ${HDF5_LIBRARIES} numgrid-objects OpenMP::OpenMP_CXX Eigen3::Eigen h5pp::h5pp)
endif()

if (BUILD_OPENCAP)
message("Building OpenCAP executable.")
find_package(Python ${PY_VERSION} REQUIRED)
add_executable(opencap ${SOURCES})
add_library(opencap-shared ${SOURCES})
target_link_libraries(opencap ${HDF5_LIBRARIES} numgrid-objects pybind11::module pybind11::embed Eigen3::Eigen h5pp::h5pp OpenMP::OpenMP_CXX)
target_link_libraries(opencap-shared ${HDF5_LIBRARIES} numgrid-objects pybind11::module pybind11::embed Eigen3::Eigen h5pp::h5pp OpenMP::OpenMP_CXX)
install(TARGETS opencap
    DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
    )
    
option(ENABLE_TESTS "Enable unit tests" ON)
if(ENABLE_TESTS)
enable_testing()
add_subdirectory(tests)
endif()

endif()
