#ifndef LEGENDREPOLY_H
#define LEGENDREPOLY_H
#include <iostream>
#include <cmath>
#include <vector>
#include <iostream>
#include <limits>

const double epsilon = std::numeric_limits<double>::epsilon();
using namespace std;

class LegendrePoly {
private:
    int maxDegree;
    int dim;
    int numBasis;
    vector<vector<double>> mesh;
    vector<vector<double>> cellsize;
    vector<vector<double>> coeff;
    void setNumBasis(const int dim, const int maxDegree);
    void setMesh(vector<vector<double>>& mesh);
public:
    long long int getNumBasis(const int dim, const int maxDegree);
    LegendrePoly() {}
    LegendrePoly(int dim, int maxDegree, vector<vector<double>>& mesh) {
        setNumBasis(dim, maxDegree);
        setMesh(mesh);
    }
    vector<vector<double>> l2Projection(double (*f)(vector<double>));
    double solve(const vector<double> x);
};
#endif