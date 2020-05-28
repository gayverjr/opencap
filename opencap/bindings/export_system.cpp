#include <pybind11/pybind11.h>
#include "System.h"

namespace py = pybind11;


PYBIND11_MODULE(pycap, m) {
    py::class_<System>(m, "System")
		.def(py::init<>())
        .def("set_geometry", &System::set_geometry);
}

