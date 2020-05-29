#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "System.h"

namespace py = pybind11;


PYBIND11_MODULE(pycap, m) {
    py::class_<System>(m, "System")
		.def(py::init<>())
        .def("set_geometry", &System::set_geometry)
		.def("print_dict",&System::print_dict)
		.def("print_dict2",&System::print_dict2)
		;
}

