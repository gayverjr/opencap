#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "System.h"
#include "ProjectedCAP.h"

namespace py = pybind11;


PYBIND11_MODULE(pycap, m) {
    py::class_<System>(m, "System")
		.def(py::init<py::dict>())
		.def("get_overlap_mat",&System::get_overlap_mat)
		;
    py::class_<Projected_CAP>(m, "Projected_CAP")
		.def(py::init<System,py::dict,size_t,std::string>())
    	.def("run",&Projected_CAP::run)
    	.def("get_ao_cap",&Projected_CAP::get_ao_cap)
    	.def("get_projected_cap",&Projected_CAP::get_projected_cap)
    	.def("compute_ao_cap",&Projected_CAP::compute_ao_cap)
    	.def("compute_projected_cap",&Projected_CAP::compute_projected_cap)
		.def("get_H",&Projected_CAP::get_H)
		.def("add_tdm",(void (Projected_CAP::*)(py::array_t<double>&,size_t,size_t))
				&Projected_CAP::add_tdm)
		.def("add_tdm",(void (Projected_CAP::*)(py::array_t<double>&,py::array_t<double>&,
				size_t,size_t)) &Projected_CAP::add_tdm)
		.def("read_data",&Projected_CAP::read_electronic_structure_data)
	;
}

