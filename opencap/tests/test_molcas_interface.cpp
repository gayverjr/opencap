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
#include <limits>
#include "Atom.h"
#include "BasisSet.h"
#include "molcas_interface.h"
#include "overlap.h"
#include "System.h"
#include "transforms.h"

TEST(RASSI_SUITE, RASSI_GEOMETRY) {
    std::vector<Atom> atms = read_geometry_from_rassi("../tests/openmolcas/cart.rassi.h5");
    ASSERT_EQ (atms.size(),2);
    ASSERT_EQ (atms[0].Z,7);
    ASSERT_EQ (atms[0].coords[0],0.00);
}

TEST(RASSI_SUITE,RASSI_CART_BASIS)
{
    std::vector<Atom> atms = read_geometry_from_rassi("../tests/openmolcas/cart.rassi.h5");
    BasisSet bs = read_basis_from_rassi("../tests/openmolcas/cart.rassi.h5",atms);
    ASSERT_EQ (bs.Nshells,16);
    ASSERT_EQ (bs.Nbasis,38);
}

TEST(RASSI_SUITE,RASSI_BASIS)
{
    std::vector<Atom> atms = read_geometry_from_rassi("../tests/openmolcas/xms_symm.rassi.h5");
    BasisSet bs = read_basis_from_rassi("../tests/openmolcas/xms_symm.rassi.h5",atms);
    ASSERT_EQ (bs.Nshells,12);
    ASSERT_EQ (bs.Nbasis,28);
}

TEST(RASSI_SUITE,RASSI_CART_OVERLAP)
{
    std::vector<Atom> atms = read_geometry_from_rassi("../tests/openmolcas/cart.rassi.h5");
    BasisSet bs = read_basis_from_rassi("../tests/openmolcas/cart.rassi.h5",atms);
    auto smat = read_rassi_overlap("../tests/openmolcas/cart.rassi.h5",bs);
    ASSERT_EQ (smat.cols(),38);
}

TEST(RASSI_SUITE,RASSI_OVERLAP)
{
    std::vector<Atom> atms = read_geometry_from_rassi("../tests/openmolcas/xms_symm.rassi.h5");
    BasisSet bs = read_basis_from_rassi("../tests/openmolcas/xms_symm.rassi.h5",atms);
    auto smat = read_rassi_overlap("../tests/openmolcas/xms_symm.rassi.h5",bs);
    auto smat2 = read_rassi_overlap("../tests/openmolcas/xms_nosymm.rassi.h5",bs);
	Eigen::MatrixXd cart_smat(bs.num_carts(),bs.num_carts());
	compute_analytical_overlap(bs,cart_smat);
	uniform_cart_norm(cart_smat,bs);
	Eigen::MatrixXd smat3(bs.Nbasis,bs.Nbasis);
	cart2spherical(cart_smat,smat3,bs);
    ASSERT_NEAR ((smat-smat2).sum(),0.0,0.0000001);
    ASSERT_NEAR ((smat-smat3).sum(),0.0,0.0000001);
}

TEST(RASSI_SUITE, RASSI_DMS)
{
    //rassi
    std::vector<Atom> atms = read_geometry_from_rassi("../tests/openmolcas/cart.rassi.h5");
    BasisSet bs = read_basis_from_rassi("../tests/openmolcas/cart.rassi.h5",atms);
    size_t nstates = 4;
    std::vector<std::vector<Eigen::MatrixXd>> alpha_dms;
	std::vector<std::vector<Eigen::MatrixXd>> beta_dms;
    read_rassi_tdms(alpha_dms,beta_dms,"../tests/openmolcas/cart.rassi.h5",bs,nstates);
    ASSERT_EQ (beta_dms[0][0].cols(),38);
    //mpssi
    std::vector<Atom> atms2 = read_geometry_from_rassi("../tests/openmolcas/dmrg_nevpt2.mpssi.h5");
    BasisSet bs2 = read_basis_from_rassi("../tests/openmolcas/dmrg_nevpt2.mpssi.h5",atms);
    std::vector<std::vector<Eigen::MatrixXd>> alpha_dms2;
	std::vector<std::vector<Eigen::MatrixXd>> beta_dms2;
    read_rassi_tdms(alpha_dms2,beta_dms2,"../tests/openmolcas/dmrg_nevpt2.mpssi.h5",bs2,nstates);
    ASSERT_EQ (beta_dms2[0][0].cols(),28);
}

TEST(MOLCAS_OUTPUT, PARSE_NEVPT2)
{
    Eigen::MatrixXd pc_nevpt2 = read_nevpt2_heff(4, "../tests/openmolcas/dmrg_nevpt2.out", "pc-nevpt2");
    Eigen::MatrixXd sc_nevpt2 = read_nevpt2_heff(4, "../tests/openmolcas/dmrg_nevpt2.out", "sc-nevpt2");
    ASSERT_EQ (sc_nevpt2(0,0), -109.2552448931);
    ASSERT_EQ (pc_nevpt2(0,0), -109.2554753551);
}

TEST(MOLCAS_OUTPUT, PARSE_CASPT2)
{
    Eigen::MatrixXd xmscaspt2 = read_mscaspt2_heff(4, "../tests/openmolcas/xms_nosymm.out"); 
    Eigen::MatrixXd mscaspt2 = read_mscaspt2_heff(4, "../tests/openmolcas/cart.out");
    Eigen::MatrixXd no_shift = read_mscaspt2_heff(5, "../tests/openmolcas/no_shift.out");
    ASSERT_EQ (mscaspt2(0,0), -109.29238639);
    ASSERT_EQ (no_shift(0,0), -95.50812900);
    //rotated matrix
    ASSERT_NEAR (xmscaspt2(3,3), -108.858019,0.00001);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}




