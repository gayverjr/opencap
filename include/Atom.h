#include <vector>
#include <array>
#include <string>
#pragma once
using namespace std;

struct element
{
        unsigned short Z;
        std::string name;
        std::string symbol;
};

class Atom {
public:
	std::string symbol;
	std::string name;
	int Z;
	std::array<double,3> coords;

public:
	void ang_to_bohr();
	Atom(string symbol,double x,double y, double z);
	element get_element_info(std::string symbol);
};
