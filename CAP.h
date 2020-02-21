/*
 * CAP.h
 *
 *  Created on: Feb 19, 2020
 *      Author: JG
 */

#ifndef CAP_H_
#define CAP_H_

#include <iostream>
#include <math.h>
#include "Atom.h"

using namespace std;

// Base class
class cap_pot {
   public:
      // pure virtual function providing interface framework.
      virtual double eval_pot(double x, double y, double z) = 0;
      virtual ~cap_pot() {}
};

class box_cap: public cap_pot {
   public:
	  double cap_x;
	  double cap_y;
	  double cap_z;
	  box_cap(double x,double y, double z):cap_x(x),cap_y(y),cap_z(z){};
      double eval_pot(double x, double y, double z)
      {
         double result = 0;
         if(abs(x)>cap_x)
        	 result += (abs(x)-cap_x) * (abs(x)-cap_x);
         if(abs(y)>cap_y)
        	 result += (abs(y)-cap_y) * (abs(y)-cap_y);
         if(abs(z)>cap_z)
        	 result += (abs(z)-cap_z) * (abs(z)-cap_z);
         return result;
      }
};

class voronoi_cap: public cap_pot {
   public:
	  double r_cut;
	  std::vector<Atom> molecule;
	  voronoi_cap(double cutoff_radius, std::vector<Atom> atoms): r_cut(cutoff_radius),molecule(atoms) {};
      double eval_pot(double x, double y, double z)
      {
          double atom_distances[molecule.size()];
          double r_closest=1000.0;
          //find r closest and fill up our distances array
          for(int j=0;j<molecule.size();j++)
          {
              if(molecule[j].Z!=0)
              {
				  double dist_x= (x-molecule[j].coords[0]) * (x-molecule[j].coords[0]);
				  double dist_y= (y-molecule[j].coords[1]) * (y-molecule[j].coords[1]);
				  double dist_z= (z-molecule[j].coords[2]) * (z-molecule[j].coords[2]);
				  double dist=sqrt(dist_x+dist_y+dist_z);
				  if(dist<r_closest)
					 r_closest=dist;
				  atom_distances[j]=dist;
              }
          }
          double weights[molecule.size()];
          for(int j=0;j<molecule.size();j++)
          {
        	  if(molecule[j].Z!=0)
        		  weights[j]=1/pow((pow(atom_distances[j],2.0)-pow(r_closest,2.0) +1),2.0);
        	  else
        		  weights[j]=0;
          }
          double numerator=0.0;
          double denominator=0.0;
          for(int j=0;j<molecule.size();j++)
          {
        	  numerator+=atom_distances[j]*atom_distances[j]*weights[j];
        	  denominator+=weights[j];
          }
          double r=sqrt(numerator/denominator);
          if(r<r_cut)
        	  return 0;
          else
        	  return (r-r_cut)*(r-r_cut);
      }
};



#endif /* CAP_H_ */
