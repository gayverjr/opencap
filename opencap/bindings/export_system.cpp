#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "System.h"
#include "ProjectedCAP.h"
#include <pybind11/numpy.h>

namespace py = pybind11;


PYBIND11_MODULE(pycap, m) {
    py::class_<System>(m, "System")
		.def(py::init<py::dict>())
		.def("get_overlap_mat",&System::get_overlap_mat)
		;
    py::class_<Projected_CAP>(m, "Projected_CAP")
		.def(py::init<System,size_t,std::string>())
    	.def("run",&Projected_CAP::run)
    	.def("get_AO_CAP",&Projected_CAP::get_AO_CAP)
		.def("get_H",&Projected_CAP::get_H)
    	.def("get_CAP_mat",&Projected_CAP::get_CAP_mat)
		.def("add_tdm",(void (Projected_CAP::*)(py::array_t<double>&,size_t,size_t))
				&Projected_CAP::add_tdm)
		.def("add_tdm",(void (Projected_CAP::*)(py::array_t<double>&,py::array_t<double>&,
				size_t,size_t)) &Projected_CAP::add_tdm)
		.def("set_h0",&Projected_CAP::set_h0)
		.def("read_data",&Projected_CAP::read_electronic_structure_data)
		.def("set_cap_params",&Projected_CAP::set_cap_params)
	;
}

