#include <iostream>
#include "./include/gaussLegendre.h"
#include <gsl/gsl_integration.h>
#include "./include/SolFunc2D.h"
double f(double x, double y) {
    return x * y;
}

int main() {
    int n = 100; // 分割数，用于数值积分
    double xa = 0;
    double xb = 2 * M_PI;
    double ya = 0;
    double yb = 2 * M_PI;
    size_t Nx = 40;
    size_t Ny = 40;
    SolFunc2D F(2, xa, xb, ya, yb, Nx, Ny);
    F.setCoeff(f);
    cout << F({ M_PI, M_PI }) << endl;
    cout << F({ 0, 0 }) << endl;
    return 0;
}