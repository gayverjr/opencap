/*
 * readMolcasHDF5.cpp
 *
 *  Created on: Mar 15, 2020
 *      Author: JG
 */
#include <armadillo>

std::vector<std::vector<arma::mat>> readRassiHDF5(std::string dmat_filename)
{
	/*
	arma::mat C = B.slice(0);
	std::cout << "Number of columns:" << C.n_cols << " Number of rows:" << sqrt(C.n_rows);
	arma::mat D = C.col(0);
	D.reshape(sqrt(C.n_rows),sqrt(C.n_rows));
	D.print();*/
	std::vector<std::vector<arma::mat>> opdms;
	//first lets load in the hdf5 file
	arma::cube rassi_data;
	rassi_data.load(arma::hdf5_name(dmat_filename, "SFS_TRANSITION_DENSITIES"));
	//now lets loop over the matrices
	for (size_t i=0;i<rassi_data.n_slices;i++)
	{
		std::vector<arma::mat> state_row;
		for(size_t j=0;j<rassi_data.n_slices;j++)
		{
			arma::mat opdm = rassi_data.slice(0).col(j);
			opdm.reshape(sqrt(opdm.n_rows),sqrt(opdm.n_rows));
			state_row.push_back(opdm);
		}
		opdms.push_back(state_row);
	}
	//for some reason, i-->j where i>j are set to 0... set them to j-->i
    for (size_t i=0;i<rassi_data.n_slices;i++)
    {
    	for(size_t j=0;j<i;j++)
    	{
    		opdms[i][j]= opdms[j][i];
    	}
    }
    return opdms;
}
