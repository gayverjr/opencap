#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "System.h"
#include "ProjectedCAP.h"
#include <pybind11/eigen.h>
#include <Eigen/Dense>

namespace py = pybind11;


PYBIND11_MODULE(pycap, m) {
    py::class_<System>(m, "System")
		.def(py::init<py::dict>())
		.def("get_overlap_mat",&System::get_overlap_mat)
		.def("check_overlap_mat",&System::check_overlap_mat,py::arg("smat"),
				py::arg("ordering"),py::arg("basis_file") = "")
		;
    py::class_<Projected_CAP>(m, "Projected_CAP")
		.def(py::init<System,py::dict,size_t,std::string>())
    	.def("run",&Projected_CAP::run)
    	.def("get_ao_cap",&Projected_CAP::get_ao_cap)
    	.def("get_projected_cap",&Projected_CAP::get_projected_cap)
    	.def("compute_ao_cap",&Projected_CAP::compute_ao_cap)
    	.def("compute_projected_cap",&Projected_CAP::compute_projected_cap)
		.def("get_H",&Projected_CAP::get_H)
		.def("add_tdm",  &Projected_CAP::add_tdm,py::arg("tdm"),
				py::arg("row_idx"),py::arg("col_idx"),
				py::arg("ordering"),py::arg("basis_file") = "")
		.def("add_tdms", &Projected_CAP::add_tdms,py::arg("alpha_density"),
				py::arg("beta_density"),py::arg("row_idx"),py::arg("col_idx"),
				py::arg("ordering"),py::arg("basis_file") = "")
		.def("read_data",&Projected_CAP::read_electronic_structure_data)
	;
}

