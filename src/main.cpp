#include <iostream>
#include "./include/gaussLegendre.h"



double f(double x, double y) {
    return x * y; // Example function
}

int main() {
    // Example usage
    double xa = 0;
    double xb = 1;
    double ya = 0;
    double yb = 1;
    int n = 5;

    auto result = gaussLegendre2D(f, { xa, xb }, { ya, yb }, n);
    std::cout << "Integral result: " << result << std::endl;

    return 0;
}