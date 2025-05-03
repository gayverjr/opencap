/*Copyright (c) 2021 James Gayvert

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

/*
 * CAP.cpp
 */

#include "CAP.h"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <sstream>
#include <string>
#include "AOCAP.h"
#include "Atom.h"
#include "BasisSet.h"
#include "gto_ordering.h"
#include "molcas_interface.h"
#include "opencap_exception.h"
#include "overlap.h"
#include "qchem_interface.h"
#include "System.h"
#include "transforms.h"
#include "utils.h"
#include "cap_types.h"


CAP::CAP(System &my_sys,std::map<std::string, std::string> params)
{
	system = my_sys;
	python = false;
	verify_method(params);
	parameters=params;
	cap_integrator = AOCAP(system.atoms,params);
	stringstream ss(parameters["nstates"]);
	ss >> nstates;
	read_in_zero_order_H();
	read_in_dms();
}

void CAP::define_cap_function(py::dict dict,const std::function<std::vector<double>(std::vector<double> &, std::vector<double> &, 
std::vector<double> &, std::vector<double> &)> &cap_func)
{
	std::vector<std::string> valid_keywords = {"cap_type","cap_x","cap_y","cap_z", "thresh",
			"r_cut","radial_precision","angular_points","do_numerical"};
	std::map<std::string, std::string> params;
    for (auto item : dict)
    {
    	std::string key = py::str(item.first).cast<std::string>();
		std::string value = py::str(item.second).cast<std::string>();
		transform(key.begin(),key.end(),key.begin(),::tolower);
		if (std::find(valid_keywords.begin(),
				valid_keywords.end(),key)==valid_keywords.end())
			opencap_throw("Invalid key in dictionary:`" + key + "'\n");
		params[key]=value;
    }
	for (auto item: params)
		parameters[item.first]=item.second;
	if(cap_func == nullptr)
		cap_integrator = AOCAP(system.atoms,params);
	else
		cap_integrator = AOCAP(system.atoms,params,cap_func);
}

CAP::CAP(System my_sys, py::dict dict, size_t num_states,const std::function<std::vector<double>(std::vector<double> &, std::vector<double> &, 
std::vector<double> &, std::vector<double> &)> &cap_func)
{
	python = true;
	system = my_sys;
	nstates = num_states;
	if(num_states<1)
		opencap_throw("Error: not enough states to run calculation.");
	define_cap_function(dict,cap_func);
}

void CAP::read_in_zero_order_H()
{
	std::string method = parameters["method"];
	transform(method.begin(),method.end(),method.begin(),::tolower);
	if (parameters.find("h0_file")!=parameters.end())
	{
		std::cout << "Reading in Hamiltonian from:" << parameters["h0_file"] << std::endl;
		ZERO_ORDER_H = read_matrix(nstates,parameters["h0_file"]);
		std::string message = "Successfully read in zeroth order Hamiltonian from file:" + parameters["h0_file"];
		if(python)
			py::print(message);
		else
			std::cout << message << std::endl;
	}
	else if (compare_strings(parameters["package"],"qchem") && parameters.find("qchem_output")!=parameters.end())
	{
		if(compare_strings(method,"eom"))
			ZERO_ORDER_H = read_qchem_eom_energies(nstates,parameters["qchem_output"]);
		else if (compare_strings(method,"tddft"))
			ZERO_ORDER_H = read_qchem_tddft_energies(nstates,parameters["qchem_output"]);
		else
			opencap_throw("Error: unsupported Q-Chem method. Currently EOM and TDDFT are supported.");
		std::string message = "Successfully read in zeroth order Hamiltonian from file:" + parameters["qchem_output"];
		if(python)
			py::print(message);
		else
			std::cout << message << std::endl;
	}
	else if (compare_strings(parameters["package"],"openmolcas") && parameters.find("molcas_output")!=parameters.end())
	{
		if(compare_strings(method,"ms-caspt2") || compare_strings(method,"xms-caspt2"))
			ZERO_ORDER_H = read_mscaspt2_heff(nstates,parameters["molcas_output"],rotation_matrix);
		else if(compare_strings(method,"sa-casscf"))
			ZERO_ORDER_H = read_sacasscf_h(nstates,parameters["molcas_output"],rotation_matrix);
		else if(compare_strings(method,"sc-nevpt2") || compare_strings(method,"pc-nevpt2"))
			ZERO_ORDER_H = read_nevpt2_heff(nstates,parameters["molcas_output"],method);
		else
			opencap_throw("Error: unsupported OpenMolcas method. Currently (X)-MS-CASPT2, and PC/SC-NEVPT2 are supported.");
		std::string message = "Successfully read in zeroth order Hamiltonian from file:" + parameters["molcas_output"];
		if(python)
			py::print(message);
		else
			std::cout << message << std::endl;
	}
	else
	{
		std::string message = "No zeroth order Hamiltonian specified. Substituting zero matrix instead.";
		if(python)
			py::print(message);
		else
			std::cout << message << std::endl;
		ZERO_ORDER_H = Eigen::MatrixXd::Zero(nstates,nstates);
	}
}

void CAP::read_in_dms()
{
	if (compare_strings(parameters["package"],"qchem"))
	{
		try
		{
			std::string message= "Reading densities from file:" + parameters["qchem_fchk"];
			if(python)
				py::print(message);
			else
				std::cout << message << std::endl;
			qchem_read_dms(alpha_dms,beta_dms,parameters["qchem_fchk"],system.bs,nstates);
			message= "Done.";
			if(python)
				py::print(message);
			else
				std::cout << message << std::endl;
		}
		catch (exception &e)
		{
			opencap_rethrow("Failed to read TDMs from Q-Chem .fchk.");
		}
	}
	else if (compare_strings(parameters["package"],"openmolcas"))
	{
		try
		{
			std::string message = "Reading densities from file:" + parameters["rassi_h5"];
			if(python)
				py::print(message);
			else
				std::cout << message << std::endl;
			read_rassi_tdms(alpha_dms,beta_dms,parameters["rassi_h5"],system.bs,nstates);
			message = "Done.";
			if(python)
				py::print(message);
			else
				std::cout << message << std::endl;
		}
		catch (exception &e)
		{
			opencap_rethrow("Failed to read TDMs from file:"+parameters["rassi_h5"]);
		}
	}
	else
		opencap_throw("Only OpenMolcas rassi.h5 and Q-Chem .fchk are supported.");

}

void CAP::compute_projected_cap()
{
	if (AO_CAP_MAT.cols()==0)
		integrate_cap();
	verify_data();
	Eigen::MatrixXd CAP_matrix(nstates,nstates);
	CAP_matrix = Eigen::MatrixXd::Zero(nstates,nstates);
	for (size_t row_idx=0;row_idx<CAP_matrix.rows();row_idx++)
	{
		for (size_t col_idx=0;col_idx<CAP_matrix.cols();col_idx++)
		{
			Eigen::MatrixXd alpha_mat_prod = alpha_dms[row_idx][col_idx]*AO_CAP_MAT;
			Eigen::MatrixXd beta_mat_prod = beta_dms[row_idx][col_idx]*AO_CAP_MAT;
			CAP_matrix(row_idx,col_idx) =  alpha_mat_prod.trace()+beta_mat_prod.trace();
			CAP_matrix(row_idx,col_idx) = -1.0* CAP_matrix(row_idx,col_idx);
		}
	}
	if (rotation_matrix.cols()!=0)
	{
		std::cout << "Warning: rotating CAP matrix (U^dagger*W*U) using the following rotation matrix U:" << std::endl << rotation_matrix << std::endl;
		CAP_matrix = rotation_matrix.transpose()*CAP_matrix*rotation_matrix;
	}
	CAP_MAT = CAP_matrix;
}

void CAP::compute_projected_capG()
{
    if (AO_CAPG_MAT.empty()) {
        integrate_capG();
    }
    verify_data();
    for (std::map<char, Eigen::MatrixXd>& matX : AO_CAPG_MAT) {
        std::map<char, Eigen::MatrixXd> CAPG_matrix;
        CAPG_matrix['x'] = Eigen::MatrixXd::Zero(nstates, nstates);
        CAPG_matrix['y'] = Eigen::MatrixXd::Zero(nstates, nstates);
        CAPG_matrix['z'] = Eigen::MatrixXd::Zero(nstates, nstates);

        for (char axis : {'x', 'y', 'z'}) {
            for (size_t row_idx = 0; row_idx < CAPG_matrix[axis].rows(); row_idx++) 
			{
                for (size_t col_idx = 0; col_idx < CAPG_matrix[axis].cols(); col_idx++) 
				{
                    Eigen::MatrixXd alpha_mat_prod = alpha_dms[row_idx][col_idx] * matX[axis];
                    Eigen::MatrixXd beta_mat_prod = beta_dms[row_idx][col_idx] * matX[axis];
                    CAPG_matrix[axis](row_idx, col_idx) = alpha_mat_prod.trace() + beta_mat_prod.trace();
                    CAPG_matrix[axis](row_idx, col_idx) = -1.0 * CAPG_matrix[axis](row_idx, col_idx);
                }
            }
            if (rotation_matrix.cols() != 0) {
                //std::cout << "Warning: rotating CAP matrix (U^dagger*W*U) using the following rotation matrix U:" << std::endl << rotation_matrix << std::endl;
                CAPG_matrix[axis] = rotation_matrix.transpose() * CAPG_matrix[axis] * rotation_matrix;
            }
        }
        CAPG_MAT.push_back(CAPG_matrix);
    }
}

void CAP::compute_projected_cap_der()
{
	if (AO_CAP_DER_MAT.empty()) {
        integrate_cap_der();
    }
    verify_data();
    for (std::map<char, Eigen::MatrixXd>& matX : AO_CAP_DER_MAT) {
        std::map<char, Eigen::MatrixXd> CAP_der_matrix;
        CAP_der_matrix['x'] = Eigen::MatrixXd::Zero(nstates, nstates);
        CAP_der_matrix['y'] = Eigen::MatrixXd::Zero(nstates, nstates);
        CAP_der_matrix['z'] = Eigen::MatrixXd::Zero(nstates, nstates);

        for (char axis : {'x', 'y', 'z'}) {
            for (size_t row_idx = 0; row_idx < CAP_der_matrix[axis].rows(); row_idx++) 
			{
                for (size_t col_idx = 0; col_idx < CAP_der_matrix[axis].cols(); col_idx++) 
				{
                    Eigen::MatrixXd alpha_mat_prod = alpha_dms[row_idx][col_idx] * matX[axis];
                    Eigen::MatrixXd beta_mat_prod = beta_dms[row_idx][col_idx] * matX[axis];
                    CAP_der_matrix[axis](row_idx, col_idx) = alpha_mat_prod.trace() + beta_mat_prod.trace();
                    CAP_der_matrix[axis](row_idx, col_idx) = -1.0 * CAP_der_matrix[axis](row_idx, col_idx);
                }
            }
            if (rotation_matrix.cols() != 0) {
                //std::cout << "Warning: rotating CAP matrix (U^dagger*W*U) using the following rotation matrix U:" << std::endl << rotation_matrix << std::endl;
                CAP_der_matrix[axis] = rotation_matrix.transpose() * CAP_der_matrix[axis] * rotation_matrix;
            }
        }
        CAP_DER_MAT.push_back(CAP_der_matrix);
    }
}

void CAP::integrate_cap()
{
    Eigen::MatrixXd cap_mat(system.bs.num_carts(),system.bs.num_carts());
    cap_mat = Eigen::MatrixXd::Zero(system.bs.num_carts(),system.bs.num_carts());
    auto start = std::chrono::high_resolution_clock::now();
	cap_integrator.compute_ao_cap_mat(cap_mat,system.bs);
    auto stop = std::chrono::high_resolution_clock::now();
    auto total_time = std::chrono::duration<double>(stop-start).count();
    if(python)
    	py::print("Integration time:"+std::to_string(total_time));
    else
    	std::cout << "Integration time:" << std::to_string(total_time) << std::endl;
    uniform_cart_norm(cap_mat,system.bs);
    Eigen::MatrixXd cap_spherical(system.bs.Nbasis,system.bs.Nbasis);
    cart2spherical(cap_mat,cap_spherical,system.bs);
    AO_CAP_MAT = cap_spherical;
}

void CAP::integrate_capG()
{
    std::vector<std::map<char, Eigen::MatrixXd>> capG_mat_MD;

    auto start = std::chrono::high_resolution_clock::now();
    cap_integrator.computeG_ao_cap_mat(capG_mat_MD, system.bs);
    auto stop = std::chrono::high_resolution_clock::now();
    auto total_time = std::chrono::duration<double>(stop - start).count();
    if (python)
        py::print("Integration time:" + std::to_string(total_time));
    else
        std::cout << "Integration time:" << std::to_string(total_time) << std::endl;

    for (const std::map<char, Eigen::MatrixXd>& matX : capG_mat_MD) {
        std::map<char, Eigen::MatrixXd> capG_spherical;
        capG_spherical['x'] = Eigen::MatrixXd::Zero(system.bs.Nbasis, system.bs.Nbasis);
        capG_spherical['y'] = Eigen::MatrixXd::Zero(system.bs.Nbasis, system.bs.Nbasis);
        capG_spherical['z'] = Eigen::MatrixXd::Zero(system.bs.Nbasis, system.bs.Nbasis);

        for (char axis : {'x', 'y', 'z'}) {
            Eigen::MatrixXd matX_copy = matX.at(axis).eval();
            uniform_cart_norm(matX_copy, system.bs);
            cart2spherical(matX_copy, capG_spherical[axis], system.bs);
        }

        AO_CAPG_MAT.push_back(capG_spherical);
    }
}

void CAP::integrate_cap_der()
{
	std::vector<std::map<char, Eigen::MatrixXd>> cap_der_mat_MD;

    auto start = std::chrono::high_resolution_clock::now();
    cap_integrator.compute_ao_cap_der_mat(cap_der_mat_MD, system.bs);
    auto stop = std::chrono::high_resolution_clock::now();
    auto total_time = std::chrono::duration<double>(stop - start).count();
    if (python)
        py::print("Integration time:" + std::to_string(total_time));
    else
        std::cout << "Integration time:" << std::to_string(total_time) << std::endl;

    for (const std::map<char, Eigen::MatrixXd>& matX : cap_der_mat_MD) {
        std::map<char, Eigen::MatrixXd> cap_der_spherical;
        cap_der_spherical['x'] = Eigen::MatrixXd::Zero(system.bs.Nbasis, system.bs.Nbasis);
        cap_der_spherical['y'] = Eigen::MatrixXd::Zero(system.bs.Nbasis, system.bs.Nbasis);
        cap_der_spherical['z'] = Eigen::MatrixXd::Zero(system.bs.Nbasis, system.bs.Nbasis);

        for (char axis : {'x', 'y', 'z'}) {
            Eigen::MatrixXd matX_copy = matX.at(axis).eval();  // Make a copy of the matrix.
            uniform_cart_norm(matX_copy, system.bs);
            cart2spherical(matX_copy, cap_der_spherical[axis], system.bs);
        }
        AO_CAP_DER_MAT.push_back(cap_der_spherical);
    }
}

void CAP::compute_cap_on_grid(py::array_t<double>& x, py::array_t<double>& y, py::array_t<double>& z, py::array_t<double>& grid_w)
{
    py::buffer_info bufx = x.request();
    py::buffer_info bufy = y.request();
    py::buffer_info bufz = z.request();
    py::buffer_info bufw = grid_w.request();
    double* ptrx = (double*)bufx.ptr;
    double* ptry = (double*)bufy.ptr;
    double* ptrz = (double*)bufz.ptr;
	double* ptrw = (double*)bufw.ptr;
	size_t num_points = bufx.size;
    Eigen::MatrixXd cap_mat(system.bs.num_carts(),system.bs.num_carts());
    cap_mat = Eigen::MatrixXd::Zero(system.bs.num_carts(),system.bs.num_carts());
    auto start = std::chrono::high_resolution_clock::now();
	cap_integrator.compute_cap_on_grid(cap_mat,system.bs,ptrx,ptry,ptrz,ptrw,num_points);
    auto stop = std::chrono::high_resolution_clock::now();
    auto total_time = std::chrono::duration<double>(stop-start).count();
    if(python)
    	py::print("Integration time:"+std::to_string(total_time));
    else
    	std::cout << "Integration time:" << std::to_string(total_time) << std::endl;
    uniform_cart_norm(cap_mat,system.bs);
    Eigen::MatrixXd cap_spherical(system.bs.Nbasis,system.bs.Nbasis);
    cart2spherical(cap_mat,cap_spherical,system.bs);
	if(AO_CAP_MAT.cols()==0)
		AO_CAP_MAT = cap_spherical;
	else
		AO_CAP_MAT += cap_spherical;
}

void CAP::compute_ao_cap(py::dict dict,const std::function<std::vector<double>(std::vector<double> &, std::vector<double> &, 
std::vector<double> &, std::vector<double> &)> &cap_func)
{
    py::print("Redefining CAP parameters...\n");
	define_cap_function(dict,cap_func);
    python = true;
    integrate_cap();
}

void CAP::check_overlap_matrix()
{
	//get overlap matrix
	if (compare_strings(parameters["package"],"qchem"))
		OVERLAP_MAT = qchem_read_overlap(parameters["qchem_fchk"],system.bs);
	else if (compare_strings(parameters["package"],"openmolcas"))
		OVERLAP_MAT = read_rassi_overlap(parameters["rassi_h5"],system.bs);
	system.check_overlap_mat(OVERLAP_MAT,"molden");
}

void CAP::renormalize_cap(Eigen::MatrixXd smat, std::string ordering,
		std::string basis_file)
{
	if(system.OVERLAP_MAT.rows() != smat.rows() || system.OVERLAP_MAT.cols() != smat.cols())
		opencap_throw("Error: Dimension of supplied overlap matrix is incorrect.");
	std::vector<bf_id> ids;
	if(compare_strings(ordering,"pyscf"))
		ids = get_pyscf_ids(system.bs);
	else if(compare_strings(ordering,"openmolcas"))
	{
		if(basis_file=="")
			opencap_throw("Error: OpenMolcas ordering requires a valid HDF5 file "
					"specified with the basis_file optional argument.");
		ids = get_molcas_ids(system.bs,basis_file);
	}
	else if(compare_strings(ordering,"qchem"))
		ids = get_qchem_ids(system.bs);
	else if(compare_strings(ordering,"psi4"))
		ids = get_psi4_ids(system.bs);
	else if(compare_strings(ordering,"molden"))
		ids = system.bs.bf_ids;
	else
		opencap_throw(ordering +" ordering is not supported.");
	to_opencap_ordering(smat,system.bs,ids);
	std::vector<double> scalars;
	for (size_t i=0;i<smat.rows();i++)
		scalars.push_back(sqrt(smat(i,i)));
	for (size_t i=0;i<AO_CAP_MAT.rows();i++)
	{
		for(size_t j=0;j<AO_CAP_MAT.cols();j++)
			AO_CAP_MAT(i,j)= AO_CAP_MAT(i,j)* scalars[i] * scalars[j];
	}
}

void CAP::renormalize()
{
	if(OVERLAP_MAT.cols()==0 && python)
		py::print("Error: no overlap matrix to use for re-normalization. Use \"read_data\" first, or"
				" use the \"renormalize_cap\" function instead.");
	else
		renormalize_cap(OVERLAP_MAT,"molden");
}

void CAP::verify_method(std::map<std::string,std::string> params)
{
	if (!python)
	{
		if(params.find("nstates")==params.end())
			opencap_throw("Error: missing the 'nstates' keyword. "
					"Please specify the number of electronic states.")
	}
	if(params.find("package")==params.end())
		opencap_throw("Error: missing the 'package' keyword. "
				"Please choose a supported package.");
	std::string package_name = params["package"];
	if(params.find("method")==params.end())
		opencap_throw("Error: missing the 'method' keyword. "
				"Please choose a supported method.");
	std::string method = params["method"];
	transform(method.begin(),method.end(),method.begin(),::tolower);
	if (compare_strings(package_name,"qchem"))
	{
		std::vector<std::string> supported = {"eom","tddft"};
		if (std::find(supported.begin(), supported.end(), method) == supported.end())
		{
			std::string message = "Error: unsupported Q-Chem method. OpenCAP currently supports: ";
			for (auto mthd: supported)
				message += mthd + ", ";
			opencap_throw(message);
		}
		if (params.find("qchem_fchk")==params.end())
			opencap_throw("Error: missing keyword: qchem_fchk.");
	}
	else if(compare_strings(package_name,"openmolcas"))
	{
		std::vector<std::string> supported = {"ms-caspt2","xms-caspt2","pc-nevpt2","sc-nevpt2", "sa-casscf"};
		if (std::find(supported.begin(), supported.end(), method) == supported.end())
		{
			std::string message = "Error: unsupported OpenMolcas method. OpenCAP currently supports: ";
			for (auto mthd: supported)
				message += mthd + ", ";
			opencap_throw(message);
		}
		if (params.find("rassi_h5")==params.end())
			opencap_throw("Error: missing keyword: rassi_h5.");
	}
	else
		opencap_throw("Error: unsupported package. Only QChem and OpenMolcas are currently supported.");
}

void CAP::run()
{
	try
	{
		if (parameters.find("ignore_overlap")==parameters.end()||parameters["ignore_overlap"]=="false")
			check_overlap_matrix();
		integrate_cap();
        renormalize();
		compute_projected_cap();
	}
	catch(exception &e)
	{
		opencap_rethrow("Error initializing CAP calculation.");
	}
}

void CAP::verify_data()
{
	//check that number of states is correct
	if(!(nstates == alpha_dms.size() && nstates == beta_dms.size()))
		opencap_throw("Error: number of states found:"+std::to_string(alpha_dms.size()) + " does not match the number of states "
				"previously specified:"+std::to_string(nstates));
	//check that dimensionality of density matrices is correct
	//alpha
	for(size_t i=0;i<nstates;i++)
	{
		for(size_t j=0;j<nstates;j++)
			if(!(alpha_dms[i][j].cols() == alpha_dms[i][j].rows() && alpha_dms[i][j].rows() == system.bs.Nbasis))
				opencap_throw("Error: Dimensionality of the density matrices do not match"
						" the basis set specified in the system object.")
	}
	//beta
	for(size_t i=0;i<nstates;i++)
	{
		for(size_t j=0;j<nstates;j++)
		{
			if(!(beta_dms[i][j].cols() == beta_dms[i][j].rows() && beta_dms[i][j].rows() == system.bs.Nbasis))
				opencap_throw("Error: Dimensionality of the density matrices do not match"
			"the basis set specified in the system object.")
		}
	}
}

Eigen::MatrixXd CAP::get_ao_cap(std::string ordering,std::string basis_file)
{
    if (AO_CAP_MAT.cols()==0)
		//std::cout << "LOOK HERE"<< AO_CAP_MAT << std::endl;
        integrate_cap();
    if(ordering!="")
    {
        Eigen::MatrixXd reordered_cap = AO_CAP_MAT;
        std::vector<bf_id> ids;
        if(compare_strings(ordering,"pyscf"))
        ids = get_pyscf_ids(system.bs);
        else if(compare_strings(ordering,"openmolcas"))
        {
            if(basis_file=="")
            opencap_throw("Error: OpenMolcas ordering requires a valid HDF5 file "
                          "specified with the basis_file optional argument.");
            ids = get_molcas_ids(system.bs,basis_file);
        }
        else if(compare_strings(ordering,"qchem"))
            ids = get_qchem_ids(system.bs);
        else if(compare_strings(ordering,"psi4"))
            ids = get_psi4_ids(system.bs);
        else if(compare_strings(ordering,"bagel"))
            ids = get_bagel_ids(system.bs);
        else if(compare_strings(ordering,"molden"))
            ids = system.bs.bf_ids;
        else
            opencap_throw("Error: " + ordering +" is unsupported.");
        reorder_matrix(reordered_cap,system.bs.bf_ids,ids);
        return reordered_cap;
    }
    return AO_CAP_MAT;
    
}
std::vector<std::map<char, Eigen::MatrixXd>> CAP::get_ao_capG(std::string ordering,std::string basis_file)
{
    if (AO_CAPG_MAT.empty()){
       	integrate_capG();
	}
	size_t matX_idx=0;
	if(ordering!="") {
		for (std::map<char, Eigen::MatrixXd>& reordered_cap : AO_CAPG_MAT){
			for (char axis : {'x', 'y', 'z'}){
        		std::vector<bf_id> ids;
        		if(compare_strings(ordering,"pyscf"))
        			ids = get_pyscf_ids(system.bs);
        		else if(compare_strings(ordering,"openmolcas"))
        		{
            		if(basis_file=="")
            		opencap_throw("Error: OpenMolcas ordering requires a valid HDF5 file "
                          	"specified with the basis_file optional argument.");
            		ids = get_molcas_ids(system.bs,basis_file);
        		}
        		else if(compare_strings(ordering,"qchem"))
            		ids = get_qchem_ids(system.bs);
        		else if(compare_strings(ordering,"psi4"))
            		ids = get_psi4_ids(system.bs);
        		else if(compare_strings(ordering,"bagel"))
            		ids = get_bagel_ids(system.bs);
        		else if(compare_strings(ordering,"molden"))
            		ids = system.bs.bf_ids;
        		else
            		opencap_throw("Error: " + ordering +" is unsupported.");
        		reorder_matrix(reordered_cap[axis],system.bs.bf_ids,ids);
			}
			AO_CAPG_MAT[matX_idx] = reordered_cap;
			matX_idx +=1;
    	}
		return AO_CAPG_MAT;
	}
	
    return AO_CAPG_MAT;
    
}
Eigen::MatrixXd CAP::get_projected_cap()
{
	return CAP_MAT;
}
std::vector<std::map<char, Eigen::MatrixXd>> CAP::get_projected_capG()
{
	return CAPG_MAT;
}
std::vector<std::map<char, Eigen::MatrixXd>> CAP::get_projected_cap_der()
{
	return CAP_DER_MAT;
}
Eigen::MatrixXd CAP::get_density(size_t row_idx, size_t col_idx, bool beta)
{
	if(alpha_dms.size()!=nstates)
		opencap_throw("Error: need to populate densities first.")
	if(row_idx >= nstates || col_idx >=nstates)
		opencap_throw("Error: invalid state index.")
	if(beta)
		return beta_dms[row_idx][col_idx];
	else
		return alpha_dms[row_idx][col_idx];
}

Eigen::MatrixXd CAP::get_H()
{
	return ZERO_ORDER_H;
}

void CAP::add_tdms(Eigen::MatrixXd &alpha_density,
		Eigen::MatrixXd &beta_density,size_t row_idx, size_t col_idx,
		std::string ordering,std::string basis_file)
{
	// if this is the first one, we need to create the data structures
	if(alpha_dms.size()!=nstates)
	{
		alpha_dms = std::vector<std::vector<Eigen::MatrixXd>>(nstates,
				std::vector<Eigen::MatrixXd>(nstates));
		beta_dms = std::vector<std::vector<Eigen::MatrixXd>>(nstates,
				std::vector<Eigen::MatrixXd>(nstates));
	}
	Eigen::MatrixXd alpha_dm = alpha_density;
	Eigen::MatrixXd beta_dm = beta_density;
	if(!(alpha_dm.cols() == alpha_dm.rows() && alpha_dm.cols() == system.bs.Nbasis &&
			beta_dm.cols() == beta_dm.rows() && beta_dm.cols() == system.bs.Nbasis))
		opencap_throw("Error: Dimensionality of the density matrices "
				"do not match the basis set specified in the system object.");
	std::vector<bf_id> ids;
	if(compare_strings(ordering,"pyscf"))
		ids = get_pyscf_ids(system.bs);
	else if(compare_strings(ordering,"openmolcas"))
	{
		if(basis_file=="")
			opencap_throw("Error: OpenMolcas ordering requires a valid HDF5 file "
					"specified with the basis_file optional argument.");
		ids = get_molcas_ids(system.bs,basis_file);
	}
	else if(compare_strings(ordering,"qchem"))
		ids = get_qchem_ids(system.bs);
	else if(compare_strings(ordering,"psi4"))
		ids = get_psi4_ids(system.bs);
	else if(compare_strings(ordering,"molden"))
		ids = system.bs.bf_ids;
	else if(compare_strings(ordering,"bagel"))
		ids = get_bagel_ids(system.bs);
	else
		opencap_throw("Error: " + ordering +" is unsupported.");
	to_opencap_ordering(alpha_dm,system.bs,ids);
	to_opencap_ordering(beta_dm,system.bs,ids);
	alpha_dms [row_idx][col_idx] = alpha_dm;
	beta_dms [row_idx][col_idx] = beta_dm;
}

Eigen::MatrixXd CAP::qcsoftware_to_opencap_ordering(Eigen::MatrixXd mat_in, std::string ordering, std::string basis_file)
{
    std::cout << "Returns " << ordering << " to OpenCAP ordering for a matrix.\n" << std::endl;
    std::vector<bf_id> input_ids;

    std::vector<std::tuple<int, int>> swap_indices;

    if (compare_strings(ordering, "pyscf"))
        input_ids = get_pyscf_ids(system.bs);
    else if (compare_strings(ordering, "openmolcas"))
    {
        if (basis_file == "")
            opencap_throw("Error: OpenMolcas ordering requires a valid HDF5 file "
                          "specified with the basis_file optional argument.");
        input_ids = get_molcas_ids(system.bs, basis_file);
    }
    else if (compare_strings(ordering, "qchem"))
        input_ids = get_qchem_ids(system.bs);
    else if (compare_strings(ordering, "psi4"))
        input_ids = get_psi4_ids(system.bs);
    else if (compare_strings(ordering, "molden"))
        input_ids = system.bs.bf_ids;
    else if (compare_strings(ordering, "bagel"))
        input_ids = get_bagel_ids(system.bs);
    else
        opencap_throw("Error: " + ordering + " is unsupported.");

    for (size_t i = 0; i < input_ids.size(); i++)
    {
        for (size_t j = 0; j < system.bs.bf_ids.size(); j++)
        {
            if (input_ids[i] == system.bs.bf_ids[j])
            {
                swap_indices.push_back(std::make_tuple(i, j));
            }
        }
    }

    Eigen::MatrixXd per_mat(system.bs.Nbasis, system.bs.Nbasis);
    per_mat = Eigen::MatrixXd::Zero(system.bs.Nbasis, system.bs.Nbasis);

    for (auto t : swap_indices)
        per_mat(std::get<0>(t), std::get<1>(t)) = 1;

    // permute indices: P^T * A * P
    Eigen::MatrixXd mat_out(system.bs.Nbasis, system.bs.Nbasis);
    mat_out = per_mat.transpose() * mat_in * per_mat;

    return mat_out;
}


void CAP::add_tdm(Eigen::MatrixXd tdm,size_t row_idx, size_t col_idx,std::string ordering,
		std::string basis_file)
{
	// if this is the first one, we need to create the data structures
	if(alpha_dms.size()!=nstates)
	{
		alpha_dms = std::vector<std::vector<Eigen::MatrixXd>>(nstates,
				std::vector<Eigen::MatrixXd>(nstates));
		beta_dms = std::vector<std::vector<Eigen::MatrixXd>>(nstates,
				std::vector<Eigen::MatrixXd>(nstates));
	}
	Eigen::MatrixXd dmat;
	dmat = 0.5 * tdm;
	if(!(dmat.cols() == dmat.rows()&& dmat.cols() == system.bs.Nbasis))
		opencap_throw("Error: Dimensionality of the density matrix "
				"does not match the basis set specified in the system object.");
	std::vector<bf_id> ids;
	if(compare_strings(ordering,"pyscf"))
		ids = get_pyscf_ids(system.bs);
	else if(compare_strings(ordering,"openmolcas"))
	{
		if(basis_file=="")
			opencap_throw("Error: OpenMolcas ordering requires a valid HDF5 file "
					"specified with the basis_file optional argument.");
		ids = get_molcas_ids(system.bs,basis_file);
	}
	else if(compare_strings(ordering,"qchem"))
		ids = get_qchem_ids(system.bs);
	else if(compare_strings(ordering,"psi4"))
		ids = get_psi4_ids(system.bs);
	else if(compare_strings(ordering,"molden"))
		ids = system.bs.bf_ids;
	else if(compare_strings(ordering,"bagel"))
		ids = get_bagel_ids(system.bs);
	else
		opencap_throw("Error: " + ordering +" is unsupported.");
	to_opencap_ordering(dmat,system.bs,ids);
	alpha_dms[row_idx][col_idx] = dmat;
	beta_dms[row_idx][col_idx] = dmat;
}



void CAP::read_electronic_structure_data(py::dict dict)
{
	std::vector<std::string> valid_keywords = {"method","qchem_output","rassi_h5","qchem_fchk","molcas_output","package","h0_file"};
	std::map<std::string, std::string> params;
    for (auto item : dict)
    {
    	std::string key = py::str(item.first).cast<std::string>();
    	std::string value = py::str(item.second).cast<std::string>();
		transform(key.begin(),key.end(),key.begin(),::tolower);
    	if (std::find(valid_keywords.begin(),
    			valid_keywords.end(),key)==valid_keywords.end())
    		opencap_throw("Invalid key in dictionary:`" + key + "'\n");
		params[key]=value;
    }
	try
	{
		verify_method(params);
		for (auto item: params)
			parameters[item.first]=item.second;
		read_in_zero_order_H();
		read_in_dms();
		check_overlap_matrix();
	}
	catch(exception &e)
	{
		for (auto item: dict)
		{
			std::string key = py::str(item.first).cast<std::string>();
	    	std::string value = py::str(item.second).cast<std::string>();
			if(parameters[key] == value)
			{
				auto it1 = parameters.find(key);
				parameters.erase(it1);
			}
		}
		throw ;
	}
}