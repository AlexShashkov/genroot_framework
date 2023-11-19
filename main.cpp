#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>

#include "framework.h"

using std::fma;

int main() {
    int l = 6;
    std::vector<double> roots(l);
    std::vector<double> a(l+1, 0.0);
    generate_polynomial(l, 0, 0, 0, 1.0, 0.0, 20.0, roots, a);

std::cout << "\nFRAMEWORK ROOTS:\n";
    for(auto &el : roots){
        std::cout << "(x-" << el << ")";
    }

std::cout << "\nFRAMEWORK COEFFS:\n";
    for(auto &el : a){
        std::cout<< el << ',';
    }   

    return 0;
}