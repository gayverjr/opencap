/*
 * CAP.h
 *
 *  Created on: Feb 19, 2020
 *      Author: JG
 */

#ifndef CAP_H_
#define CAP_H_

#include <iostream>
#include "Atom.h"

using namespace std;

// Base class
class cap_func {
   public:
      // pure virtual function providing interface framework.
      virtual double eval_pot(double x, double y, double z) = 0;
      virtual ~cap_func() {}
};

class box_cap: public cap_func {
   public:
	  box_cap(double cap_x,double cap_y, double cap_z);
      double eval_pot(double x, double y, double z) {
         return 1.0;
      }
};

class voronoi_cap: public cap_func {
   public:
	  voronoi_cap(double r_cut, std::vector<Atom> atoms);
      double eval_pot(double x, double y, double z) {
         return 1.0;
      }
};



#endif /* CAP_H_ */
