#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "System.h"
#include "ProjectedCAP.h"

namespace py = pybind11;


PYBIND11_MODULE(pycap, m) {
    py::class_<System>(m, "System")
		.def(py::init<py::dict>());
    py::class_<Projected_CAP>(m, "Projected_CAP")
		.def(py::init<System>())
    	.def("run",&Projected_CAP::run);
}

