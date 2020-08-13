#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "System.h"
#include "CAP.h"
#include <pybind11/eigen.h>
#include <Eigen/Dense>

namespace py = pybind11;


PYBIND11_MODULE(pyopencap_cpp, m) {
    py::class_<System>(m, "System")
		.def(py::init<py::dict>(),"Constructs System object.")
		.def("get_overlap_mat",&System::get_overlap_mat,"Returns overlap matrix.")
		.def("check_overlap_mat",&System::check_overlap_mat,py::arg("smat"),
				py::arg("ordering"),py::arg("basis_file") = "","Compares input overlap matrix to "
						"internal overlap to check basis set ordering.")
		.def("get_basis_ids",&System::get_basis_ids,"Returns a string of the basis function ids. "
				"Each ID has the following format:"
				"atom index,shell number,l,m")
		;
    py::class_<CAP>(m, "CAP")
		.def(py::init<System,py::dict,size_t,std::string>(),"Constructs CAP object.")
    	.def("get_ao_cap",&CAP::get_ao_cap, "Returns CAP matrix in AO basis.")
    	.def("get_perturb_cap",&CAP::get_perturb_cap, "Returns CAP matrix in state basis.")
    	.def("compute_ao_cap",&CAP::compute_ao_cap, "Computes CAP matrix in AO basis.")
    	.def("compute_perturb_cap",&CAP::compute_perturb_cap, "Computes CAP matrix in state basis using"
    			" transition density matrices.")
		.def("get_H",&CAP::get_H, "Returns zeroth order Hamiltonian read from file.")
		.def("add_tdm",  &CAP::add_tdm,py::arg("tdm"),
				py::arg("initial_idx"),py::arg("final_idx"),
				py::arg("ordering"),py::arg("basis_file") = "",
				"Adds spin-traced tdm to CAP object at specified indices. The optional argument basis_file"
				" is required when using the OpenMolcas interface, and it must point to the path to the rassi.5 file.")
		.def("add_tdms", &CAP::add_tdms,py::arg("alpha_density"),
				py::arg("beta_density"),py::arg("initial_idx"),py::arg("final_idx"),
				py::arg("ordering"),py::arg("basis_file") = "",
				"Adds alpha/beta tdms to CAP object at specified indices. The optional argument basis_file"
				" is required when using the OpenMolcas interface, and it must point to the path to the rassi.5 file.")
		.def("read_data",&CAP::read_electronic_structure_data, py::arg("es_dict"), "Reads electronic structure data "
				"specified in dictionary.")
		.def("renormalize_cap",&CAP::renormalize_cap, py::arg("smat"),
				py::arg("ordering"),py::arg("basis_file") = "","Re-normalizes AO CAP matrix using input overlap"
						"matrix.")
		.def("renormalize",&CAP::renormalize,"Re-normalizes AO CAP using electronic structure data.")
	;
}

