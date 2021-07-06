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
#include <omp.h>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>

#include "Atom.h"
#include "BasisSet.h"
#include "CAP.h"
#include "gto_ordering.h"
#include "InputParser.h"
#include "keywords.h"
#include "molden_parser.h"
#include "System.h"
#include "utils.h"
#include "TrajectoryAnalysis.h"

TEST(CAP, TEST_VORONOI)
{
   std::tuple<System,std::map<std::string,std::string>> inp_data = parse_input("../tests/cap/test_voronoi.in");
	std::map<std::string,std::string> params = std::get<1>(inp_data); 
   CAP pc(std::get<0>(inp_data),get_params_for_field(params,"projected_cap"));
   pc.integrate_cap();
   params["r_cut"] = "4.00";
   CAP pc2(std::get<0>(inp_data),get_params_for_field(params,"projected_cap"));
   pc2.integrate_cap();
   ASSERT_NE(pc2.AO_CAP_MAT.sum(),pc.AO_CAP_MAT.sum());
}

TEST(CAP, TEST_BOX)
{
   std::tuple<System,std::map<std::string,std::string>> inp_data = parse_input("../tests/cap/test_box.in");
	std::map<std::string,std::string> params = std::get<1>(inp_data); 
   CAP pc(std::get<0>(inp_data),get_params_for_field(params,"projected_cap"));
   pc.integrate_cap();
   params["cap_z"] = "5.00";
   CAP pc2(std::get<0>(inp_data),get_params_for_field(params,"projected_cap"));
   pc2.integrate_cap();
   ASSERT_NE(pc2.AO_CAP_MAT.sum(),pc.AO_CAP_MAT.sum());
}

TEST(CAP, TEST_GRID)
{
   std::tuple<System,std::map<std::string,std::string>> inp_data = parse_input("../tests/cap/test_box.in");
	std::map<std::string,std::string> params = std::get<1>(inp_data); 
   CAP pc(std::get<0>(inp_data),get_params_for_field(params,"projected_cap"));
   pc.integrate_cap();
   params["radial_precision"] = "12";
   CAP pc2(std::get<0>(inp_data),get_params_for_field(params,"projected_cap"));
   pc2.integrate_cap();
   ASSERT_NE(pc2.AO_CAP_MAT.sum(),pc.AO_CAP_MAT.sum());
   params["angular_points"] = "350";
   CAP pc3(std::get<0>(inp_data),get_params_for_field(params,"projected_cap"));
   pc3.integrate_cap();
   ASSERT_NE(pc2.AO_CAP_MAT.sum(),pc3.AO_CAP_MAT.sum());
}

TEST(CAP, TEST_PARALLEL)
{
   std::tuple<System,std::map<std::string,std::string>> inp_data = parse_input("../tests/cap/test_box.in");
	std::map<std::string,std::string> params = std::get<1>(inp_data); 
   CAP pc(std::get<0>(inp_data),get_params_for_field(params,"projected_cap"));
   pc.run();
   omp_set_num_threads(1);
   CAP pc2(std::get<0>(inp_data),get_params_for_field(params,"projected_cap"));
   pc2.run();
   ASSERT_NEAR(pc2.AO_CAP_MAT.sum(),pc.AO_CAP_MAT.sum(),0.000001);
   ASSERT_NEAR(pc2.CAP_MAT.sum(),pc.CAP_MAT.sum(),0.000001);
}

TEST(CAP,TEST_TRAJ)
{
   std::tuple<System,std::map<std::string,std::string>> inp_data = parse_input("../tests/cap/test_box.in");
	std::map<std::string,std::string> params = std::get<1>(inp_data); 
   CAP pc(std::get<0>(inp_data),get_params_for_field(params,"projected_cap"));
   pc.run();
   CAPHamiltonian caph(pc.ZERO_ORDER_H,pc.CAP_MAT,get_params_for_field(params,"trajectory"),"../tests/cap/test_box.in");
	caph.run_trajectory();
   SUCCEED();
}

TEST(CAP,BAD_DIFFUSE)
{
    std::tuple<System,std::map<std::string,std::string>> inp_data = parse_input("../tests/cap/test_bad_diffuse.in");
    std::map<std::string,std::string> params = std::get<1>(inp_data);
    CAP pc(std::get<0>(inp_data),get_params_for_field(params,"projected_cap"));
    try
    {
        pc.integrate_cap();
    }
    catch (const std::exception&)
    {
        SUCCEED();
    }
}

TEST(CAP,TEST_IMPORTH)
{
    auto h0_full = read_matrix(10,"../tests/cap/heff.in");
    ASSERT_EQ(h0_full(0,0),-109.312105);
    auto h0_diag = read_matrix(5,"../tests/cap/heff_diag.in");
    ASSERT_EQ(h0_diag(0,0),-109.313539);
}

int main(int argc, char **argv) 
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
