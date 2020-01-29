#include <iostream>
#include <vector>
#include <math.h>
#include "Atom.h"
using namespace std;

Atom::Atom(string sym,double x_coord,double y_coord, double z_coord)
{
	symbol=sym;
	coords.push_back(x_coord);
	coords.push_back(y_coord);
	coords.push_back(z_coord);
}
