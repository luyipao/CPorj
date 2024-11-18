#include <iostream>
#include <gsl/gsl_integration.h>
#include "./include/SolFunc2D.h"
double f(double x, double y) {
    return sin(x) * sin(y);
}
vector<vector<VectorXd>> L(const vector<vector<VectorXd>>& coeff, const MatrixXd& B, const MatrixXd& D, const MatrixXd& E) {
    int N_x = coeff.size() - 1;
    int N_y = coeff[0].size() - 1;

    vector<vector<VectorXd>> L(N_x + 1, vector<VectorXd>(N_y + 1, VectorXd(coeff[0][0].size())));
    for (int i = 1; i <= N_x; ++i) {
        for (int j = 1; j <= N_y; ++j) {
            VectorXd temp = -B * coeff[i][j] - D * coeff[i - 1][j] - E * coeff[i][j - 1];
            L[i][j] = temp;
        }
    }

    // periodic boundary condition
    for (int j = 1;j < N_y + 1; ++j) {
        L[0][j] = L[N_x][j];
    }
    for (int i = 0; i < N_x + 1; ++i) {
        L[i][0] = L[i][N_y];
    }
    return L;
}

vector<vector<VectorXd>> RK(const vector<vector<VectorXd>>& coeff, const MatrixXd& B, const MatrixXd& D, const MatrixXd& E, double t) {
    int N = coeff.size() - 1;
    int basisNum = coeff[0][0].size();
    // 计算 k1
    auto k1 = L(coeff, B, D, E);
    for (int i = 0; i < N + 1; ++i) {
        for (int j = 0; j < N + 1; ++j) {
            k1[i][j] = coeff[i][j] + t * k1[i][j];
        }
    }
    cout << "k1[0][1] : \n" << k1[0][1] << endl;
    cout << "k1[1][0] : \n" << k1[1][0] << endl;
    cout << "k1[1][1] : \n" << k1[1][1] << endl;
    auto k2 = L(k1, B, D, E);
    for (int i = 0; i < N + 1; ++i) {
        for (int j = 0; j < N + 1; ++j) {
            k2[i][j] = 3.0 / 4.0 * coeff[i][j] + 1.0 / 4.0 * k1[i][j] + 1.0 / 4 * t * k2[i][j];
        }
    }
    // vector<vector<VectorXd>> k3(N, vector<VectorXd>(N, VectorXd(basisNum)));
    auto k3 = L(k2, B, D, E);
    for (int i = 0; i < N + 1; ++i) {
        for (int j = 0; j < N + 1; ++j) {
            k3[i][j] = 1.0 / 3.0 * coeff[i][j] + 2.0 / 3.0 * k2[i][j] + 2.0 / 3.0 * t * k3[i][j];
        }
    }
    return k3;
}

int main() {

    double CFL = 0.05;
    double xa = 0;
    double xb = 2 * M_PI;
    double ya = 0;
    double yb = 2 * M_PI;
    double tb = 4 * M_PI;
    int Nx = 40;
    int Ny = 40;
    int k = 2;
    double dx = (xb - xa) / Nx;
    double dy = (yb - ya) / Ny;

    SolFunc2D P(2, xa, xb, ya, yb, Nx, Ny);
    cout << "-----------计算系数-----------" << endl;
    P.setCoeff(f);
    cout << endl;

    for (int i = 0; i <= Nx; ++i) {
        cout << P({ i * dx, M_PI / 2 }) << ' ';
    }
    cout << endl;

    int coeffSize = getBasisNum(2, k);

    // construct mass matrix
    MatrixXd B(coeffSize, coeffSize);
    MatrixXd D(coeffSize, coeffSize);
    MatrixXd E(coeffSize, coeffSize);
    for (int i = 0; i < coeffSize; ++i) {
        for (int j = 0; j < coeffSize; ++j) {
            auto m = getVarDegree(2, i);
            auto n = getVarDegree(2, j);

            // 使用 () 访问矩阵元素
            B(i, j) = sqrt((2 * m[1] + 1) * (2 * n[1] + 1) / dx / dy) * (m[0] == n[0]);
            B(i, j) += sqrt((2 * m[0] + 1) * (2 * n[0] + 1)) / dx * (m[1] == n[1]);
            B(i, j) -= (m[0] > n[0]) * ((m[0] - n[0]) % 2 == 1) * 2 * (m[1] == n[1]) * sqrt((2 * m[0] + 1) * (2 * n[0] + 1)) / dx;
            B(i, j) -= (m[1] > n[1]) * ((m[1] - n[1]) % 2 == 1) * 2 * (m[0] == n[0]) * sqrt((2 * m[1] + 1) * (2 * n[1] + 1)) / dy;

            // 使用 std::pow 进行幂运算
            D(i, j) = -std::pow(-1, m[0]) * (m[1] == n[1]) * sqrt((2 * m[0] + 1) * (2 * n[0] + 1)) / dx;
            E(i, j) = -std::pow(-1, m[1]) * (m[0] == n[0]) * sqrt((2 * m[1] + 1) * (2 * n[1] + 1)) / dy;
        }
    }
    cout << "B:\n" << B << std::endl;

    std::cout << "D:\n" << D << std::endl;
    cout << "E:\n" << E << std::endl;

    double t = CFL * pow(sqrt(dx * dy), (k + 1) / 3.0);
    vector<double> T;
    T.emplace_back(0);

    // 得到初值函数的系数
    vector<vector<VectorXd>> coeff(Nx + 1, vector<VectorXd>(Ny + 1, VectorXd(coeffSize)));
    for (int i = 1; i < Nx + 1; ++i) {
        for (int j = 1; j < Ny + 1; ++j) {
            coeff[i][j] = P.getCoeff()[i - 1][j - 1];
        }
    }
    for (int j = 1;j < Ny + 1; ++j) {
        coeff[0][j] = coeff[Nx][j];
    }
    for (int i = 0; i < Nx + 1; ++i) {
        coeff[i][0] = coeff[i][Ny];
    }

    cout << "coeff[0][0] : " << coeff[0][0] << endl;

    cout << "-----------RK迭代-----------" << endl;
    while (T.back() < tb) {
        T.emplace_back(T.back() + t);
        coeff = RK(coeff, B, D, E, t);
        showProgressBar(tb, T.back());
    }
    cout << endl;

    vector<vector<VectorXd>> C(Nx, vector<VectorXd>(Ny, VectorXd(coeffSize)));
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            C[i][j] = coeff[i + 1][j + 1];
        }
    }
    P.setCoeff(C);
    cout << P({ 0, 0 }) << endl;
    cout << P({ M_PI,M_PI }) << endl;
    return 0;
}