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
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>

#include "Atom.h"
#include "BasisSet.h"
#include "gto_ordering.h"
#include "InputParser.h"
#include "molden_parser.h"
#include "System.h"
#include "utils.h"


TEST(GTO_ORDERING, TEST_PSI4_SPH)
{
   //spherical
	std::tuple<System,std::map<std::string,std::string>> inp_data = parse_input("../tests/system/sph_inline.in");
   System my_sys = std::get<0>(inp_data);
   auto smat = read_matrix(my_sys.bs.Nbasis,"../tests/system/psi4_sph_smat.mat");
   my_sys.check_overlap_mat(smat,"psi4");
   SUCCEED();  
}

TEST(GTO_ORDERING, TEST_PSI4_CART)
{
   std::tuple<System,std::map<std::string,std::string>> inp_data = parse_input("../tests/system/cart_inline.in");
   System my_sys = std::get<0>(inp_data);
   auto smat = read_matrix(my_sys.bs.Nbasis,"../tests/system/psi4_cart_smat.mat");
   my_sys.check_overlap_mat(smat,"psi4");
   SUCCEED();  
}

TEST(GTO_ORDERING, TEST_PYSCF_SPH)
{
   //spherical
	std::tuple<System,std::map<std::string,std::string>> inp_data = parse_input("../tests/system/test_pyscf_sph.in");
   System my_sys = std::get<0>(inp_data);
   auto smat = read_matrix(my_sys.bs.Nbasis,"../tests/system/pyscf_sph_smat.mat");
   my_sys.check_overlap_mat(smat,"pyscf");
   SUCCEED();  
}

TEST(GTO_ORDERING, TEST_PYSCF_CART)
{
   std::tuple<System,std::map<std::string,std::string>> inp_data = parse_input("../tests/system/test_pyscf_cart.in");
   System my_sys = std::get<0>(inp_data);
   auto smat = read_matrix(my_sys.bs.Nbasis,"../tests/system/pyscf_cart_smat.mat");
   my_sys.check_overlap_mat(smat,"pyscf");
   SUCCEED();  
}

TEST(MOLECULE_SPEC, TEST_INLINE_AU)
{
	std::tuple<System,std::map<std::string,std::string>> inp_data = parse_input("../tests/system/sph_inline.in");
   System my_sys = std::get<0>(inp_data);
   ASSERT_NEAR(my_sys.atoms[0].coords[2],1.037,0.0001);
}

TEST(MOLECULE_SPEC, TEST_INLINE_ANG)
{
	std::tuple<System,std::map<std::string,std::string>> inp_data = parse_input("../tests/system/cart_inline.in");
   System my_sys = std::get<0>(inp_data);
   ASSERT_NEAR(my_sys.atoms[0].coords[2],1.037,0.0001);
}

TEST(MOLECULE_SPEC, CART_BF)
{
	std::tuple<System,std::map<std::string,std::string>> inp_data = parse_input("../tests/system/cart_inline.in");
   System my_sys = std::get<0>(inp_data);
   ASSERT_EQ(my_sys.bs.Nbasis,53);
}

int main(int argc, char **argv) 
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
