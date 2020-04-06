#include <vector>
#include <array>
#pragma once
using namespace std;

class Atom {
public:
	std::string element;
	int Z;
	std::array<double,3> coords;

public:
	void ang_to_bohr();
	Atom(string symbol,double x,double y, double z);
};
