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
add_subdirectory(asa)
SET(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

find_package(OpenMP REQUIRED)
SET(H5PP_ENABLE_EIGEN3 ON)

FetchContent_Declare(
    h5pp
    GIT_REPOSITORY https://github.com/DavidAce/h5pp.git
    GIT_TAG 0d1cb5ebb6805b3d00fb323a31e91fb5998c7853
)
FetchContent_MakeAvailable(h5pp)

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG v1.14.0
)
FetchContent_MakeAvailable(googletest)


FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11
    GIT_TAG v2.13.6
)
FetchContent_MakeAvailable(pybind11)

# Have to use old populate pattern for now until we fork numgrid 
# or find another grid provider
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
target_link_libraries(pyopencap_cpp PRIVATE ${HDF5_LIBRARIES} numgrid-objects OpenMP::OpenMP_CXX Eigen3::Eigen h5pp::h5pp asa239)
endif()

if (BUILD_OPENCAP)
message("Building OpenCAP executable.")
find_package(Python ${PY_VERSION} REQUIRED)
add_executable(opencap ${SOURCES})
add_library(opencap-shared ${SOURCES})
target_link_libraries(opencap ${HDF5_LIBRARIES} numgrid-objects pybind11::module pybind11::embed Eigen3::Eigen h5pp::h5pp OpenMP::OpenMP_CXX asa239)
target_link_libraries(opencap-shared ${HDF5_LIBRARIES} numgrid-objects pybind11::module pybind11::embed Eigen3::Eigen h5pp::h5pp OpenMP::OpenMP_CXX asa239)
install(TARGETS opencap
    DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
    )
    
option(ENABLE_TESTS "Enable unit tests" ON)
if(ENABLE_TESTS)
enable_testing()
add_subdirectory(tests)
endif()

endif()
