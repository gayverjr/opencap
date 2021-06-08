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

#include "Atom.h"
#include "BasisSet.h"
#include "molden_parser.h"


TEST(MOLDEN_PARSER, PARSE_GEOMETRY_ANG)
{
   std::vector<Atom> atms = read_geometry_from_molden("../tests/molden/H2.molden");
   ASSERT_EQ (atms.size(),4);
   ASSERT_EQ (atms[0].Z,0);
   ASSERT_NEAR(atms[0].coords[2],-1.398397332,0.00001);
}

TEST(MOLDEN_PARSER, PARSE_GEOMETRY_AU)
{
   std::vector<Atom> atms = read_geometry_from_molden("../tests/molden/N2.molden");
   ASSERT_EQ (atms.size(),2);
   ASSERT_EQ (atms[0].Z,7);
   ASSERT_EQ(atms[0].coords[2],1.03934937);
}

TEST(MOLDEN_PARSER, PARSE_BASIS_CART)
{
   std::vector<Atom> atms = read_geometry_from_molden("../tests/molden/H2.molden");
   BasisSet bs = read_basis_from_molden("../tests/molden/H2.molden",atms);
   ASSERT_EQ (bs.Nshells,44);
   ASSERT_EQ (bs.Nbasis,124);
}

TEST(MOLDEN_PARSER, PARSE_BASIS_SPH)
{
   std::vector<Atom> atms = read_geometry_from_molden("../tests/molden/N2.molden");
   BasisSet bs = read_basis_from_molden("../tests/molden/N2.molden",atms);
   ASSERT_EQ (bs.Nshells,12);
   ASSERT_EQ (bs.Nbasis,28);
}

int main(int argc, char **argv) 
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
