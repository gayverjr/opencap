file(COPY qchem DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY openmolcas DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
message(STATUS "building tests")

file(GLOB files "test_*.cpp")

foreach(file ${files})
	string(REGEX REPLACE "(^.*/|\\.[^.]*$)" "" file_without_ext ${file})
	add_executable(${file_without_ext} ${file})
	target_link_libraries(${file_without_ext} opencap-shared numgrid-objects pybind11::module pybind11::embed Eigen3::Eigen h5pp::h5pp)
	add_test(${file_without_ext} ${file_without_ext})
endforeach()

