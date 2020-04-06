#include <iostream>
#include <vector>
#include <math.h>
#include "Atom.h"
#include "elements.h"



Atom::Atom(string sym,double x_coord,double y_coord, double z_coord)
{
	element=sym;
	coords = {{x_coord,y_coord,z_coord}};
	Z = 7;
}

void Atom::ang_to_bohr()
{
	const double ANG_TO_BOHR_CONV=1/0.5291772108;
	for(size_t i=0;i<coords.size();i++)
		coords[i]=coords[i]*ANG_TO_BOHR_CONV;
}
