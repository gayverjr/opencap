#include <vector>
#pragma once
using namespace std;

class Atom {
public:
	string element;
	std::vector<double> coords;
	const double ANG_TO_BOHR = 1.88973;

public:
	Atom(string symbol,double x,double y, double z);

};
