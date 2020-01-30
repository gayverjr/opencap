#include <iostream>
#include "BasisFunction.h"
#include <vector>
#include "overlap.h"
#include "BasisSet.h"
using namespace std;

int main(int argc, char **argv) {

	cout << "Hello world" << std::endl;
	//cout << overlap_integral(bf3,bf2) << std::endl;
	BasisSet bs("N2.xyz","bas.bas");
	cout << overlap_integral(bs.basis[0],bs.basis[0]) << std::endl;
	return 0;
}


