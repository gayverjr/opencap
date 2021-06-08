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

#include <gtest/gtest.h>
#include <Eigen/Dense>

#include "Atom.h"
#include "BasisSet.h"
#include "overlap.h"
#include "qchem_interface.h"
#include "System.h"
#include "transforms.h"

TEST(FCHK_SUITE, FCHK_GEOMETRY) {
    std::vector<Atom> atms = read_geometry_from_fchk("../tests/qchem/EA.fchk");
    ASSERT_EQ (atms.size(),3);
    ASSERT_EQ (atms[2].Z,0);
    ASSERT_EQ (atms[0].Z,7);
    ASSERT_EQ (atms[0].coords[0],0.00);
}

TEST(FCHK_SUITE,FCHK_CART_BASIS)
{
    std::vector<Atom> atms = read_geometry_from_fchk("../tests/qchem/EE_cart.fchk");
    BasisSet bs = read_basis_from_fchk("../tests/qchem/EE_cart.fchk",atms);
    ASSERT_EQ (bs.Nshells,15);
    ASSERT_EQ (bs.Nbasis,53);
}

TEST(FCHK_SUITE,FCHK_BASIS)
{
    std::vector<Atom> atms = read_geometry_from_fchk("../tests/qchem/EA.fchk");
    BasisSet bs = read_basis_from_fchk("../tests/qchem/EA.fchk",atms);
    ASSERT_EQ (bs.Nshells,15);
    ASSERT_EQ (bs.Nbasis,43);
}

TEST(FCHK_SUITE,FCHK_OVERLAP)
{
    std::vector<Atom> atms = read_geometry_from_fchk("../tests/qchem/EA.fchk");
    BasisSet bs = read_basis_from_fchk("../tests/qchem/EA.fchk",atms);
    auto smat = qchem_read_overlap("../tests/qchem/EA.fchk",bs);
	Eigen::MatrixXd cart_smat(bs.num_carts(),bs.num_carts());
	compute_analytical_overlap(bs,cart_smat);
	uniform_cart_norm(cart_smat,bs);
	Eigen::MatrixXd smat2(bs.Nbasis,bs.Nbasis);
	cart2spherical(cart_smat,smat2,bs);
    ASSERT_NEAR ((smat-smat2).sum(),0.0,0.0000001);
}

TEST(FCHK_SUITE, FCHK_DMS)
{
    //eomea
    std::vector<Atom> atms = read_geometry_from_fchk("../tests/qchem/EA.fchk");
    BasisSet bs = read_basis_from_fchk("../tests/qchem/EA.fchk",atms);
    size_t nstates = 3;
    std::vector<std::vector<Eigen::MatrixXd>> alpha_dms;
	std::vector<std::vector<Eigen::MatrixXd>> beta_dms;
    qchem_read_dms(alpha_dms,beta_dms,"../tests/qchem/EA.fchk",bs,nstates);
    ASSERT_EQ (beta_dms[0][0].cols(),43);
    // eomee
    read_geometry_from_fchk("../tests/qchem/EE.fchk");
    read_basis_from_fchk("../tests/qchem/EE.fchk",atms);
    std::vector<std::vector<Eigen::MatrixXd>> alpha_dms2;
	std::vector<std::vector<Eigen::MatrixXd>> beta_dms2;
    qchem_read_dms(alpha_dms2,beta_dms2,"../tests/qchem/EE.fchk",bs,nstates);
    ASSERT_EQ (beta_dms2[0][0].cols(),43);
    //tddft
    nstates = 4;
    std::vector<std::vector<Eigen::MatrixXd>> alpha_dms3;
	std::vector<std::vector<Eigen::MatrixXd>> beta_dms3;
    qchem_read_dms(alpha_dms3,beta_dms3,"../tests/qchem/tddft.fchk",bs,nstates);
    ASSERT_EQ (beta_dms3[0][0].cols(),43);
}

TEST(QCHEM_OUTPUT, PARSE_TDDFT)
{
    Eigen::MatrixXd tddft = read_qchem_tddft_energies(4,"../tests/qchem/tddft.out");
    ASSERT_EQ (tddft(0,0), -109.22753742);
}

TEST(QCHEM_OUTPUT, PARSE_EOM)
{
    Eigen::MatrixXd eomea = read_qchem_eom_energies(3, "../tests/qchem/EA.out"); 
    Eigen::MatrixXd eomee = read_qchem_eom_energies(3, "../tests/qchem/EE.out");   
    ASSERT_EQ (eomea(0,0), -108.90131733);
    ASSERT_EQ (eomee(0,0), -108.72192494); 
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}




