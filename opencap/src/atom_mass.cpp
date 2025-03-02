/*Copyright (c) 2023 Soubhik Mondal

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "atom_mass.h"
#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>

const std::unordered_map<std::string, double> MASSES = {
    {"H", 1.007825}, {"He", 4.002603}, {"Li", 7.016004}, {"Be", 9.012182},
    {"B", 11.009305}, {"C", 12.000000}, {"N", 14.003074}, {"O", 15.994915},
    {"F", 18.998403}, {"Ne", 19.992440}, {"Na", 22.989770}, {"Mg", 23.985042},
    {"Al", 26.981538}, {"Si", 27.976927}, {"P", 30.973762}, {"S", 31.972071},
    {"Cl", 34.968853}, {"Ar", 39.962383}, {"K", 38.963707}, {"Ca", 39.962591},
    {"Sc", 44.955910}, {"Ti", 47.947947}, {"V", 50.943964}, {"Cr", 51.940512},
    {"Mn", 54.938050}, {"Fe", 55.934942}, {"Co", 58.933200}, {"Ni", 57.935348},
    {"Cu", 62.929601}, {"Zn", 63.929147}, {"Ga", 68.925581}, {"Ge", 73.921178},
    {"As", 74.921596}, {"Se", 79.916522}, {"Br", 78.918338}, {"Kr", 83.911507},
    {"Rb", 84.911789}, {"Sr", 87.905614}, {"Y", 88.905848}, {"Zr", 89.904704},
    {"Nb", 92.906378}, {"Mo", 97.905408}, {"Tc", 98.907216}, {"Ru", 101.904350},
    {"Rh", 102.905504}, {"Pd", 105.903483}, {"Ag", 106.905093}, {"Cd", 113.903358},
    {"In", 114.903878}, {"Sn", 119.902197}, {"Sb", 120.903818}, {"Te", 129.906223},
    {"I", 126.904468}, {"Xe", 131.904154}, {"Cs", 132.905447}, {"Ba", 137.905241},
    {"La", 138.906348}, {"Ce", 139.905435}, {"Pr", 140.907648}, {"Nd", 141.907719},
    {"Pm", 144.912744}, {"Sm", 151.919729}, {"Eu", 152.921227}, {"Gd", 157.924101},
    {"Tb", 158.925343}, {"Dy", 163.929171}, {"Ho", 164.930319}, {"Er", 165.930290},
    {"Tm", 168.934211}, {"Yb", 173.938858}, {"Lu", 174.940768}, {"Hf", 179.946549},
    {"Ta", 180.947996}, {"W", 183.950933}, {"Re", 186.955751}, {"Os", 191.961479},
    {"Ir", 192.962924}, {"Pt", 194.964774}, {"Au", 196.966552}, {"Hg", 201.970626},
    {"Tl", 204.974412}, {"Pb", 207.976636}, {"Bi", 208.980383}, {"Po", 208.982416},
    {"At", 209.987131}, {"Rn", 222.017570}, {"Fr", 223.019731}, {"Ra", 226.025403},
    {"Ac", 227.027747}, {"Th", 232.038050}, {"Pa", 231.035879}, {"U", 238.050783}
};

// List of element symbols by atomic number (1-indexed)
const std::vector<std::string> periodic_table = {
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
};

// Map to get atomic number from symbol
const std::unordered_map<std::string, int> atomic_number_map = [] {
    std::unordered_map<std::string, int> m;
    for (size_t i = 0; i < periodic_table.size(); ++i)
        m[periodic_table[i]] = i + 1;
    return m;
}();

// Function to get element symbol from atomic number
std::string get_atom_symbol(int atomic_number) {
    if (atomic_number >= 1 && atomic_number <= static_cast<int>(periodic_table.size())) {
        return periodic_table[atomic_number - 1];
    }
    return "Invalid atomic number";
}

// Function to get atomic number from element symbol
int get_atomic_number(const std::string& atom_symbol) {
    auto it = atomic_number_map.find(atom_symbol);
    return (it != atomic_number_map.end()) ? it->second : -1;
}


