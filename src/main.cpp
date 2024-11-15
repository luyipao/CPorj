#include "include/LegendrePoly.h"
#include "include/getVarDegree.h"
#include <cmath>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

vector<vector<VectorXd>> L(vector<vector<VectorXd>>& coeff, const MatrixXd& B, const MatrixXd& D, const MatrixXd& E) {
    int N_x = coeff.size() - 1;
    int N_y = coeff[0].size() - 1;

    vector<vector<VectorXd>> L(N_x + 1, vector<VectorXd>(N_y + 1, VectorXd(coeff[0][0].size())));
    for (int j = N_y; j >= 1; --j) {
        for (int i = N_x; i >= 1; --i) {
            L[i][j] = B * coeff[i][j] + D * coeff[i - 1][j] + E * coeff[i][j - 1];
        }
    }
    // periodic boundary condition
    for (int j = 0;j < N_y + 1; ++j) {
        L[0][j] = L[N_x][j];
    }
    for (int i = 0; i < N_x + 1; ++i) {
        L[i][0] = L[i][N_y];
    }
    return L;
}

vector<vector<VectorXd>> RK(vector<vector<VectorXd>>& coeff, const MatrixXd& B, const MatrixXd& D, const MatrixXd& E, double t) {
    int N = coeff.size() - 1;
    int basisNum = coeff[0][0].size();
    // 计算 k1
    auto k1 = L(coeff, B, D, E);
    for (int i = 0; i < N + 1; ++i) {
        for (int j = 0; j < N + 1; ++j) {
            k1[i][j] = coeff[i][j] + t * k1[i][j];
        }
    }

    auto k2 = L(k1, B, D, E);
    for (int i = 0; i < N + 1; ++i) {
        for (int j = 0; j < N + 1; ++j) {
            k2[i][j] = 3.0 / 4.0 * coeff[i][j] + 1.0 / 4.0 * k1[i][j] + 1.0 / 4 * t * k2[i][j];
        }
    }
    //    vector<vector<VectorXd>> k3(N, vector<VectorXd>(N, VectorXd(basisNum)));
    auto k3 = L(k2, B, D, E);
    for (int i = 0; i < N + 1; ++i) {
        for (int j = 0; j < N + 1; ++j) {
            k3[i][j] = 1.0 / 3.0 * coeff[i][j] + 2.0 / 3.0 * k2[i][j] + 2.0 / 3.0 * t * k3[i][j];
        }
    }
    return k3;
}



int main() {
    double xa = 0;
    double xb = 2 * M_PI;
    double ya = 0;
    double yb = 2 * M_PI;
    double tb = M_PI;
    int N = 4;
    int k = 2;
    double dx = (xb - xa) / N;
    double dy = (yb - ya) / N;
    vector<double> X;
    vector<double> Y;
    for (int i = 0; i < N; ++i) {
        X.emplace_back(xa + i * dx); // 修正拼写错误
        Y.emplace_back(ya + i * dy); // 修正拼写错误
    }
    vector<vector<double>> mesh = { X, Y };
    LegendrePoly P(2, 3, mesh);
    P.getCoeff();

    int coeffDim = P.numBasis;
    // cell number
    int N_x = mesh[0].size() - 1;
    int N_y = mesh[1].size() - 1;

    vector<vector<VectorXd>> coeff(N_x + 1, vector<VectorXd>(N_y + 1, VectorXd(coeffDim)));
    for (int i = 1; i < N_x + 1; ++i) {
        for (int j = 1; j < N_y + 1; ++j) {
            coeff[i][j] = Map<VectorXd>(P.coeff[i - 1][j - 1].data(), coeffDim);
        }
    }
    // periodic boundary condition
    for (int j = 0;j < N_y + 1; ++j) {
        coeff[0][j] = coeff[N_x][j];
    }
    for (int i = 0; i < N_x + 1; ++i) {
        coeff[i][0] = coeff[i][N_y];
    }

    // construct mass matrix
    MatrixXd B(coeffDim, coeffDim);
    MatrixXd D(coeffDim, coeffDim);
    MatrixXd E(coeffDim, coeffDim);
    for (int i = 1; i < coeffDim + 1; ++i) {
        for (int j = 1; j < coeffDim + 1; ++j) {
            auto m = getVarDegree(2, i);
            auto n = getVarDegree(2, j);

            // 使用 () 访问矩阵元素
            B(i - 1, j - 1) = -sqrt((2 * m[1] + 1) * (2 * n[1] + 1) / dx / dy) * (m[0] == n[0]);
            B(i - 1, j - 1) += -(n[0] > m[0]) * ((n[0] - m[0]) % 2 == 1) * 2 * (m[1] == n[1]);
            B(i - 1, j - 1) += -(n[1] > m[1]) * ((n[1] - m[1]) % 2 == 1) * 2 * (m[0] == n[0]);

            // 使用 std::pow 进行幂运算
            D(i - 1, j - 1) = std::pow(-1, m[0] + n[0] + 1) * (m[1] == n[1]);
            E(i - 1, j - 1) = std::pow(-1, m[1] + n[1]) * (m[0] == n[0]);
        }
    }
    double CFL = 0.05;
    // 使用 std::pow 进行幂运算
    double t = CFL * (sqrt(dx * dy)) * std::pow((k + 1) / 3.0, 1.0); // 确保除法是浮点数
    vector<double> T;
    T.emplace_back(0);

    auto C = P.coeff;
    while (T.back() < tb) {
        T.emplace_back(T.back() + t);
        coeff = RK(coeff, B, D, E, t);
    }
    P.coeff = C;
    cout << P({0, 0}) << endl;
    cout << P({M_PI,M_PI}) << endl;

    return 0;
}