#include <iostream>
#include "draw.cpp"
#include "../include/gaussLegendre.h"
#include "./DDGIC4DD.cpp"
using namespace Eigen;
using namespace std;

int main() {
    int k = 1;
    int N = 100;
    double Xa = 0.0;
    double Xb = 0.6;
    double h = (Xb - Xa) / N;
    double beta0 = 1;
    double beta1 = 1.0 / 12;
    double T = 0.9;
    double CFL = 2;
    VectorXd mesh = VectorXd::LinSpaced(N + 1, Xa, Xb); // Mesh points
    MatrixXd n_C = MatrixXd::Zero(k + 1, N); // Coefficient matrix
    MatrixXd n_d_C = MatrixXd::Zero(k + 1, N); // Derivative coefficient matrix

    // construct initial function coeff
    MatrixXd coeff(k + 1, N); // 每个单元的系数向量
    for (int j = 0; j < N; ++j) {
        // 构建质量矩阵
        MatrixXd M = MatrixXd::Zero(k + 1, k + 1);
        for (int i = 0; i <= k; ++i) {
            for (int j = 0; j <= k; ++j) {
                M(i, j) = h / (i + j + 1);
            }
        }

        int gslNum = 10; // 高斯积分点个数
        auto [nodes, weights] = gaussLegendrePoints(mesh(j), mesh(j + 1), gslNum);

        VectorXd b = VectorXd::Zero(k + 1);
        for (int i = 0; i <= k; ++i) {
            b(i) = (weights.array() * n_d(nodes).array() * ((nodes.array() - mesh(j)) / h).pow(i)).sum();
        }
        // 求解线性系统
        n_C.col(j) = M.colPivHouseholderQr().solve(b);
    }
    DDGIC4DD ddgic = DDGIC4DD(N, Xa, Xb, k, beta0, beta1, T, CFL);
    ddgic.setn(n_C);
    ddgic.setn_d(n_C);
    auto start = std::chrono::high_resolution_clock::now();
    ddgic.RKDDG();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::micro> duration = end - start;
    std::cout << "耗时: " << duration.count() / 1000000.0 << " 秒\n";
    draw(ddgic.C);
}