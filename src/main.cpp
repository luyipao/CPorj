#include "include/LegendrePoly.h"
#include "include/getVarDegree.h"
#include <cmath>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

int main() {
    vector<vector<double>> mesh = { {-1, 0,1},{-1, 0,1} };
    LegendrePoly P(2, 3, mesh);
    P.getCoeff();

    int coeffDim = P.numBasis;
    // cell number
    int N_x = mesh[0].size() - 1;
    int N_y = mesh[1].size() - 1;

    vector<vector<VectorXd>> coeff(N_x + 1, vector<VectorXd>(N_y + 1, VectorXd(coeffDim)));
    for (int i = 0; i < N_x + 1; ++i) {
        for (int j = 0; j < N_y + 1; ++j) {
            coeff[i][j] = Map<VectorXd>(P.coeff[i][j].data(), coeffDim);
        }
    }


    return 0;
}