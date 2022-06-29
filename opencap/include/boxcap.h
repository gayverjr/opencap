#include <vector>
#include "Shell.h"

double boxcap(gto gto1, gto gto2, double boxlength[3],
	      int normalized = 0, int verbose = 0);
double integrate_box_cap(Shell shell1, Shell shell2, 
std::array<size_t,3> l1, std::array<size_t,3> l2,double boxlength[3],double thresh);