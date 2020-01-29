#include <vector>
#pragma once
using namespace std;

class Atom {
public:
	string symbol;
	std::vector<double> coords;

public:
	Atom(string symbol,double x,double y, double z);

};
