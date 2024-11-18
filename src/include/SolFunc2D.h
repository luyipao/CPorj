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
#include "gaussLegendre.h"
#include "showProgressBar.h"
using namespace std;
using namespace Eigen;

class SolFunc2D {
private:
    size_t maxDegree;
    size_t basisNum;

    vector<VectorXd> grid;
    vector<size_t> cellNum;
    void setNumBasis(const size_t maxDegree);
    void setGrid(const vector<double>& x, const vector<double>& y);

    vector<vector<VectorXd>> coeff;    
public:
    SolFunc2D() {}
    SolFunc2D(size_t maxDegree, double xa, double xb, double ya, double yb, size_t Nx, size_t Ny);
    SolFunc2D(int maxDegree, const vector<double>& x, const vector<double>& y) {
        setNumBasis(maxDegree);
        setGrid(x, y);
    }
    void setCoeff(function<double(double, double)> f);
    void setCoeff(const vector<vector<VectorXd>>& coeff) {
        this->coeff = coeff;
    }
    vector<vector<VectorXd>> getCoeff() {
        return this->coeff;
    }
    double operator()(const vector<double> x);

};
#endif