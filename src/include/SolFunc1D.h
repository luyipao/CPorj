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

class SolFunc1D
private:
    size_t maxDegree;
    size_t basisNum;

    VectorXd mesh;
    int cellNum;
    void setNumBasis(const size_t maxDegree);
    void setMesh(const vector<double>& x, const vector<double>& y);

    vector<VectorXd> coeff;
public:
    

#endif SOLFUNC2D_H