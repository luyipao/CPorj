#include "../include/LegendrePoly.h"
#include <iomanip>

using namespace std;
using namespace std::chrono;

double normlizedP(int n, double x, double a, double b) {
    double map_x = (2 * x - a - b) / (b - a);
    double res = sqrt((2 * n + 1) / (b - a)) * gsl_sf_legendre_Pl(n, map_x);
    return res;
}

double f(double* x) {
    return x[0] * x[1];
}

struct Params {
    vector<int> varDegree;
    vector<double> X;
    vector<double> Y;
};


double g(double* x, size_t dim, void* param) {
    Params* p = static_cast<Params*>(param);
    return f(x) *
        normlizedP(p->varDegree[0], x[0], p->X[0], p->X[1]) *
        normlizedP(p->varDegree[1], x[1], p->Y[0], p->Y[1]);
}

/**
 * @brief set information of the basis functions
 *
 * @param maxDegree basis function max degree
 * @param dim number of variables in the polynomial
 * @return int number of basis functions
 */
void LegendrePoly::setNumBasis(const int dim, const int maxDegree) {
    this->maxDegree = maxDegree;
    this->dim = dim;
    this->numBasis = getNumBasis(dim, maxDegree);
}
/**
 * @brief set the mesh,
 *
 * @param mesh v[i] represents the i-th variable's mesh, i = 0, 1, ..., dim - 1
 */
void LegendrePoly::setMesh(vector<vector<double>>& mesh) {
    for (auto v : mesh) {
        if (v.size() < 2) {
            cerr << "if your mesh have less than two points in some dimension, you should try dimensionality reduction " << endl;
        }
    }
    if (mesh.size() != this->dim) {
        cerr << "mesh size wrong" << endl;
    }

    this->mesh = mesh;
    this->cellsize = vector<vector<double>>(this->dim);
    for (int i = 0; i < this->dim; ++i) {
        for (int j = 1; j < mesh[i].size(); ++j) {
            cellsize.emplace_back(mesh[i][j] - mesh[i][j - 1]);
        }
    }
}

/**
 * @brief get the number of basis functions
 *
 * @param maxDegree basis function max degree
 * @param dim number of variables in the polynomial
 * @return int number of basis functions
 */
long long int LegendrePoly::getNumBasis(const int dim, const int maxDegree) {
    long long int a = 1;
    long long int b = 1;
    for (int i = 1; i <= maxDegree; i++) {
        a *= (dim + maxDegree);
        b *= i;
    }
    return a / b;
}


void LegendrePoly::getCoeff() {
    // only support 2d
    if (this->dim != 2) {
        cout << "dim must be 2" << endl;
    }
    // some gsl monte carlo parameters
    double res, err;
    size_t calls = 50000;
    gsl_rng_env_setup();

    const gsl_rng_type* T;
    gsl_rng* r;

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    // get the coefficients
    // initialize 
    vector<vector<vector<double>>> coeff(mesh[0].size() - 1, vector<vector<double>>(mesh[1].size() - 1, vector<double>(this->numBasis)));
    Params params;

    gsl_monte_function G = { &g, 2, &params };
    gsl_monte_plain_state* s = gsl_monte_plain_alloc(this->dim);

    for (int i = 0; i < mesh[0].size() - 1; ++i) {
        for (int j = 0; j < mesh[1].size() - 1; ++j) {
            vector<double> tempCoeff;
            for (int k = 1; k <= this->numBasis; ++k) {

                // update the parameters
                auto varDegree = getVarDegree(this->dim, k);
                params = {
                    varDegree,
                    { mesh[0][i], mesh[0][i + 1] },
                    { mesh[1][j], mesh[1][j + 1] }
                };

                double xl[this->dim] = { mesh[0][i], mesh[1][j] };
                double xu[this->dim] = { mesh[0][i + 1], mesh[1][j + 1] };


                gsl_monte_plain_integrate(&G, xl, xu, this->dim, calls, r, s, &res, &err);


                tempCoeff.emplace_back(res);
            }
            coeff[i][j] = tempCoeff;
        }
    }

    this->coeff = coeff;
    gsl_monte_plain_free(s);
    gsl_rng_free(r);
}


double LegendrePoly::operator()(vector<double> x) {

    for (int i = 0; i < this->dim; ++i) {
        if (x[i] < mesh[i].front() && x[i] >= mesh[i].back()) {
            cerr << "x isn't in the mesh" << endl;
        }
    }

    vector<int> position(dim);
    for (int i = 0; i < this->dim; ++i) {
        auto l = lower_bound(mesh[i].begin(), mesh[i].end(), x[i]);
        position[i] = l - mesh[i].begin();
        if (position[i]) {
            position[i] -= 1;
        }
    }

    double res = 0;
    for (int i = 1; i <= this->numBasis; ++i) {
        vector<int> varDegree = getVarDegree(this->dim, i);

        res += this->coeff[position[0]][position[1]][i - 1] *
            normlizedP(varDegree[0], x[0], this->mesh[0][position[0]], this->mesh[0][position[0] + 1]) *
            normlizedP(varDegree[1], x[1], this->mesh[1][position[1]], this->mesh[1][position[1] + 1]);
    }
    return res;
}