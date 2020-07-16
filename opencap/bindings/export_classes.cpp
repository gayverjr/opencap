#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "System.h"
#include "ProjectedCAP.h"
#include <pybind11/eigen.h>
#include <eigen3/Eigen/Dense>

namespace py = pybind11;


PYBIND11_MODULE(pycap, m) {
    py::class_<System>(m, "System")
		.def(py::init<py::dict>())
		.def("get_overlap_mat",&System::get_overlap_mat,"Returns the overlap matrix.")
		.def("check_overlap_mat",&System::check_overlap_mat,py::arg("smat"),
				py::arg("ordering"),py::arg("basis_file") = "","Compares input overlap matrix to "
						"internal overlap to check basis set ordering.")
		;
    py::class_<Projected_CAP>(m, "Projected_CAP")
		.def(py::init<System,py::dict,size_t,std::string>())
    	.def("run",&Projected_CAP::run, "Runs projected CAP calculation.")
    	.def("get_ao_cap",&Projected_CAP::get_ao_cap, "Returns CAP matrix in AO basis.")
    	.def("get_projected_cap",&Projected_CAP::get_projected_cap, "Returns CAP matrix in correlated basis.")
    	.def("compute_ao_cap",&Projected_CAP::compute_ao_cap, "Computes CAP matrix in AO basis.")
    	.def("compute_projected_cap",&Projected_CAP::compute_projected_cap, "Computes CAP matrix in correlated basis using"
    			"transition density matrices.")
		.def("get_H",&Projected_CAP::get_H, "Returns zeroth order Hamiltonian read from file.")
		.def("add_tdm",  &Projected_CAP::add_tdm,py::arg("tdm"),
				py::arg("row_idx"),py::arg("col_idx"),
				py::arg("ordering"),py::arg("basis_file") = "",
				"Adds spin-traced tdm to Projected CAP object at specified indices.")
		.def("add_tdms", &Projected_CAP::add_tdms,py::arg("alpha_density"),
				py::arg("beta_density"),py::arg("row_idx"),py::arg("col_idx"),
				py::arg("ordering"),py::arg("basis_file") = "",
				"Adds alpha/beta tdms to Projected CAP object at specified indices.")
		.def("read_data",&Projected_CAP::read_electronic_structure_data, py::arg("es_dict"), "Reads electronic structure data "
				"specified in dictionary.")
	;
}

