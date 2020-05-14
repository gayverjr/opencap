 /*! \file qchem_interface.h
     \brief Functions pertaining to the Q-Chem interface.
 */
#ifndef QCHEM_INTERFACE_H_
#define QCHEM_INTERFACE_H_
#include <armadillo>

/*! Read in TDMs from Q-Chem .fchk file
    \param dmat_filename: file location of the Q-Chem .fchk file
    \param nstates: number of states
    \param num_bf: number of basis functions
    \return An array of size 2, each entry holding a 2D vector corresponding to the state densities.
     The first entry is the alpha densities, the second the beta densities.
*/
std::array<std::vector<std::vector<arma::mat>>,2> qchem_read_in_dms(std::string dmat_filename,size_t nstates, size_t num_bf);
arma::mat qchem_read_overlap(std::string dmat_filename, size_t num_bf);
/*! Read in energies from Q-Chem .fchk file
    \param nstates: number of states,
    \param method: electronic structure method method,
    \param output_file: qchem .fchk file location
    \return Armadillo matrix containing the diagonal H0 matrix. Dimension is MxM, where M is number of states.
*/
arma::mat read_qchem_energies(size_t nstates,std::string method,std::string output_file);

#endif /* QCHEM_INTERFACE_H_*/
