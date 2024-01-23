// Example of how to use framework

#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>

#include "framework.h"

using std::fma;

constexpr int exponent = 2; // Exponent for big number notation
constexpr int mantissa = 4; // Mantissa for big number notation

int main() {
    int l = 5;
    std::vector<double> roots(l, 0.0);
    std::vector<double> a(l+1, 0.0);
    generate_polynomial<double, exponent, mantissa>(l, 0, l, 0,
        1e-5, -1.0, 1.0, roots, a);

std::cout << "\nFRAMEWORK ROOTS:\n";
    for(auto &el : roots){
        std::cout << "(x-(" << el << "))";
    }

std::cout << "\nFRAMEWORK COEFFS:\n";
    for(auto &el : a){
        std::cout<< el << ',';
    }   

    return 0;
}