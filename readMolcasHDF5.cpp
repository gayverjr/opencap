/*
 * readMolcasHDF5.cpp
 *
 *  Created on: Mar 15, 2020
 *      Author: JG
 */
#include <armadillo>

std::array<std::vector<std::vector<arma::mat>>,2> read_rassi_HDF5(std::string dmat_filename)
{
	std::vector<std::vector<arma::mat>> alpha_opdms;
	std::vector<std::vector<arma::mat>> beta_opdms;
	//first lets load in the hdf5 file
	arma::cube rassi_data, spin_density;
	rassi_data.load(arma::hdf5_name(dmat_filename, "SFS_TRANSITION_DENSITIES"));
	spin_density.load(arma::hdf5_name(dmat_filename, "SFS_TRANSITION_SPIN_DENSITIES"));
	//now lets loop over the matrices
	for (size_t i=0;i<rassi_data.n_slices;i++)
	{
		std::vector<arma::mat> alpha_state_row;
		std::vector<arma::mat> beta_state_row;
		for(size_t j=0;j<rassi_data.n_slices;j++)
		{
			arma::mat opdm = rassi_data.slice(i).col(j);
			arma::mat spin = spin_density.slice(i).col(i);
			spin.reshape(sqrt(spin.n_rows),sqrt(spin.n_rows));
			opdm.reshape(sqrt(opdm.n_rows),sqrt(opdm.n_rows));
			alpha_state_row.push_back(0.5*(opdm+spin));
			beta_state_row.push_back(0.5*(opdm-spin));
		}
		alpha_opdms.push_back(alpha_state_row);
		beta_opdms.push_back(beta_state_row);
	}
	//symmetric
    for (size_t i=0;i<rassi_data.n_slices;i++)
    {
    	for(size_t j=0;j<i;j++)
    	{
    		alpha_opdms[i][j]= alpha_opdms[j][i];
    		beta_opdms[i][j]= beta_opdms[j][i];
    	}
    }
    return {alpha_opdms,beta_opdms};
}
