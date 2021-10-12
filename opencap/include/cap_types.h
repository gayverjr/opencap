#include <vector>
#include <functional>
#include <map>
#include <string>
#include "Atom.h"

struct box_cap
{
    double cap_x;
    double cap_y;
    double cap_z;
    box_cap(double x, double y, double z);
    std::vector<double> operator()(std::vector<double> &x, std::vector<double> &y, std::vector<double>  &z, 
                    std::vector<double>  &grid_w);
};

// Smooth Voronoi CAP
// Thommas Sommerfeld and Masahiro Ehara
// DOI: 10.1021/acs.jctc.5b00465
struct voronoi_cap
{
    double r_cut;
    std::vector<Atom> atoms;
    voronoi_cap(double cutoff,std::vector<Atom> geometry);
    std::vector<double> operator()(std::vector<double> &x, std::vector<double> &y, std::vector<double>  &z, 
                    std::vector<double> &grid_w);
};


