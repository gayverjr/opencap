cmake_minimum_required(VERSION 3.10)
project(OpenCAP
        LANGUAGES CXX C Fortran
        )
include_directories(include)
file(GLOB SOURCES "src/*.cpp")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
include(FetchContent)

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

add_executable(OpenCAP ${SOURCES})
target_link_libraries(OpenCAP ${ARMADILLO_LIBRARIES} ${HDF5_LIBRARIES} numgrid-objects)

install(TARGETS OpenCAP
    DESTINATION ${CMAKE_INSTALL_PREFIX}
    )
