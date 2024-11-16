#ifndef SOLFUNC2D_H
#define SOLFUNC2D_H
#include <iostream>
#include <cmath>
#include <vector>
#include <iostream>
#include <chrono>
#include <algorithm>
#include <limits>
#include <functional>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_sf_legendre.h>
#include <Eigen/Dense>
#include "getVarDegree.h"
#include "getBasisNum.h"

using namespace std;
using namespace Eigen;

class SolFunc2D {
private:
    int maxDegree;
    int basisNum;

    vector<VectorXd> grid;
    vector<size_t> cellNum;

    void setNumBasis(const int maxDegree);
    void setGrid(const vector<double>& x, const vector<double>& y);
    void setCoeff(function<double(double)> f);
public:
    vector<vector<VectorXd>> coeff;
    SolFunc2D() {}
    SolFunc2D(int maxDegree, const vector<double>& x, const vector<double>& y) {
        setNumBasis(maxDegree);
        setGrid(x, y);
    }
    double operator()(vector<double> x);
};
#endif