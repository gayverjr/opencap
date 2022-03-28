add_library (asa239 asa239.cpp)
set_property(TARGET asa239 PROPERTY POSITION_INDEPENDENT_CODE ON)

# Make sure the compiler can find include files for our Hello library
# when other libraries or executables link to Hello
target_include_directories (asa239 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})