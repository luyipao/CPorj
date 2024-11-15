#ifndef LEGENDREPOLY_H
#define LEGENDREPOLY_H
#include <iostream>
#include <cmath>
#include <vector>
#include <iostream>
#include <chrono>
#include <algorithm>
#include <limits>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_sf_legendre.h>
#include "getVarDegree.h"

const double epsilon = std::numeric_limits<double>::epsilon();
using namespace std;

class LegendrePoly {
private:
    int maxDegree;
    int dim;

    vector<vector<double>> mesh;
    vector<vector<double>> cellsize;
    vector<int> cellNum;

    void setNumBasis(const int dim, const int maxDegree);
    void setMesh(vector<vector<double>>& mesh);
public:
    int numBasis;
    vector<vector<vector<double>>> coeff;
    long long int getNumBasis(const int dim, const int maxDegree);
    LegendrePoly() {}
    LegendrePoly(int dim, int maxDegree, vector<vector<double>>& mesh) {
        setNumBasis(dim, maxDegree);
        setMesh(mesh);
    }
    void getCoeff();
    // Overload operator()
    double operator()(vector<double> x);
};
#endif