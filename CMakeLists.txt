cmake_minimum_required(VERSION 3.10)
project(OpenCAP
        LANGUAGES CXX C
        )
include_directories(include)
file(GLOB SOURCES "src/*.cpp")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
find_package(armadillo)
find_package(hdf5)

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
target_link_libraries(OpenCAP armadillo hdf5 numgrid)

install(TARGETS OpenCAP
    DESTINATION ${CMAKE_INSTALL_PREFIX}
    )