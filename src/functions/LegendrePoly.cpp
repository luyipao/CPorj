#include "../include/LegendrePoly.h"

using namespace std;

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
 * @brief set the mesh, add 烧饼 in the back
 *
 * @param mesh v[i] represents the i-th variable's mesh, i = 0, 1, ..., dim, while last element represent 烧饼
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
        mesh[i].push_back(mesh[i].back() + epsilon);
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
        a *= (dim + maxDegree - i);
        b *= i;
    }
    return a / b;
}
/*
double LegendrePoly::solve(const vector<double> x) {
    if (x.size() != this->dim) {
        cerr << "point dim don't correct" << endl;
    }

    for (int i = 0; i < this->dim; ++i) {
        if (x[i] < mesh[i].first() && x[i] > mesh[i].end()) {
            cerr << "x isn't in the mesh" << endl;
        }
    }

    vector<int> position(dim);
    for (int i = 0; i < this->dim; ++i) {
        auto l = lower_bound(mesh[i].first(), mesh[i].end(), x[i]);
        if (l == mesh[i].end()) {
            cerr << "x isn't in the mesh" << endl;
        }
        position[i] = l - mesh[i].begin();
    }

    vector<vector<double>> basisValue(this->dim + 1, vector<double>(this->maxDegree + 1));
    for (int i = 0; i < this->dim; ++i) {
        for (int j = 0; j < mesh[i].size() - 1; ++i) {

        }
    }
    return 0;
}

*/

