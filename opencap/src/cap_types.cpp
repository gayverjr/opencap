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

#include "cap_types.h"
#include <math.h>
#include "opencap_exception.h"
#include <iostream>

voronoi_cap::voronoi_cap(double cutoff,std::vector<Atom> geometry)
{
    r_cut=cutoff;
    atoms=geometry;
}

std::vector<double> voronoi_cap::operator()(std::vector<double> &x, std::vector<double> &y, std::vector<double>  &z, 
                    std::vector<double>  &grid_w)
{
    std::vector<double> cap_values(x.size());
    for(size_t i=0;i<x.size();i++)
    {
        size_t num_atoms = atoms.size();
        double atom_distances[num_atoms];
        double r_closest=1000;
        //find r closest and fill up our distances array
        for(size_t j=0;j<num_atoms;j++)
        {
            if(atoms[j].Z!=0)
            {
                double dist_x= (x[i]-atoms[j].coords[0]) * (x[i]-atoms[j].coords[0]);
                double dist_y= (y[i]-atoms[j].coords[1]) * (y[i]-atoms[j].coords[1]);
                double dist_z= (z[i]-atoms[j].coords[2]) * (z[i]-atoms[j].coords[2]);
                double dist=sqrt(dist_x+dist_y+dist_z);
                if(dist<r_closest || j==0)
                    r_closest=dist;
                atom_distances[j]=dist;
            }
        }
        double weights[num_atoms];
        for(size_t q=0;q<num_atoms;q++)
        {
        if(atoms[q].Z!=0)
        {
            double weight = atom_distances[q]*atom_distances[q]-r_closest*r_closest+1;
            weights[q] = 1/(weight*weight);
        }
        else
            weights[q]=0;
        }
        double numerator=0.0;
        double denominator=0.0;
        for(size_t k=0;k<num_atoms;k++)
        {
        numerator+=atom_distances[k]*atom_distances[k]*weights[k];
        denominator+=weights[k];
        }
        double r=sqrt(numerator/denominator);
        double result=0;
        if(r>r_cut)
        result = (r-r_cut)*(r-r_cut);
        cap_values[i] = result * grid_w[i];
    }
    return cap_values;
}

box_cap::box_cap(double x,double y, double z)
{
    cap_x=x;
    cap_y=y;
    cap_z=z;
}

std::vector<double> box_cap::operator()(std::vector<double> &x, std::vector<double> &y, std::vector<double>  &z, 
                    std::vector<double>  &grid_w)
{
    std::vector<double> cap_values(x.size());
    for(size_t i=0;i<x.size();i++)
	{
		double result = 0;
		if(abs(x[i])>cap_x)
		 result += (abs(x[i])-cap_x) * (abs(x[i])-cap_x);
		if(abs(y[i])>cap_y)
		 result += (abs(y[i])-cap_y) * (abs(y[i])-cap_y);
		if(abs(z[i])>cap_z)
		 result += (abs(z[i])-cap_z) * (abs(z[i])-cap_z);
		cap_values[i] = result * grid_w[i];
	}
    return cap_values;
}