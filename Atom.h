#include <vector>
#include <array>
#pragma once
using namespace std;

class Atom {
public:
	std::string element;
	std::array<double,3> coords;
	const double ANG_TO_BOHR = 1.88973;

public:
	Atom(string symbol,double x,double y, double z);
};
