cmake_minimum_required(VERSION 3.10)
project(opencap
        LANGUAGES CXX C 
        )
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
include_directories(include)
file(GLOB SOURCES "src/*.cpp")
file(GLOB BIND_SOURCES "bindings/*.cpp")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
include(FetchContent)

FetchContent_Declare(
    carma
    GIT_REPOSITORY https://github.com/RUrlus/carma.git
)
FetchContent_GetProperties(carma)
if(NOT carma_POPULATED)
    FetchContent_Populate(carma)
    add_subdirectory(${carma_SOURCE_DIR} ${carma_BINARY_DIR})
endif()

FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11
)
FetchContent_GetProperties(pybind11)
if(NOT pybind11_POPULATED)
    FetchContent_Populate(pybind11)
    add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
endif()

FetchContent_Declare(
  numgrid
  GIT_REPOSITORY https://github.com/dftlibs/numgrid.git
)
FetchContent_GetProperties(numgrid)
if(NOT numgrid_POPULATED)
  FetchContent_Populate(numgrid)
  add_subdirectory(
    ${numgrid_SOURCE_DIR}
    ${numgrid_BINARY_DIR}
    )
endif()

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

pybind11_add_module(pycap ${SOURCES} ${BIND_SOURCES})
target_link_libraries(pycap PRIVATE ${ARMADILLO_LIBRARIES} ${HDF5_LIBRARIES} numgrid-objects carma)

add_executable(opencap ${SOURCES})
target_link_libraries(opencap ${ARMADILLO_LIBRARIES} ${HDF5_LIBRARIES} ${PYTHON_LIBRARIES} numgrid-objects pybind11::module pybind11::embed carma)

install(TARGETS opencap
    DESTINATION ${CMAKE_INSTALL_PREFIX}
    )
