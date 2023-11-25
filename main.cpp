#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>

#include "framework.h"

using std::fma;

int main() {
    int l = 5;
    std::vector<long double> roots(l, 0.0);
    std::vector<long double> a(l+1, 0.0);
    generate_polynomial(l, 0, 0, 0, 1.0L, 0.0L, 20.0L, roots, a);

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