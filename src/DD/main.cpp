#include <iostream>
#include "draw.cpp"
#include "../include/gaussLegendre.h"
#include "./DDGIC4DD.cpp"
using namespace Eigen;
using namespace std;

int main() {
    int k =2;
    double Xa = 0.0;
    double Xb = 0.6;

    double beta0 = 2;
    double beta1 = 1 / 12;
    double T = 0.3;
    double CFL = 1;



    // construct initial function coeff
    vector<int> Ns = { 20, 40, 80, 160};
    vector<VectorXd> h_values;
    vector<double> L2_errors(Ns.size(), 0.0);
    vector<double> LInf_errors(Ns.size(), 0.0);
    for (auto& N : Ns) {
        MatrixXd coeff(k + 1, N); // 每个单元的系数向量
        double h = (Xb - Xa) / N;
        VectorXd mesh = VectorXd::LinSpaced(N + 1, Xa, Xb); // Mesh points
        MatrixXd n_C = MatrixXd::Zero(k + 1, N); // Coefficient matrix
        MatrixXd n_d_C = MatrixXd::Zero(k + 1, N); // Derivative coefficient matrix
        for (int j = 0; j < N; ++j) {
            // 构建质量矩阵
            MatrixXd M = MatrixXd::Zero(k + 1, k + 1);
            for (int i = 0; i <= k; ++i) {
                for (int l = 0; l <= k; ++l) {
                    M(i, l) = h / (i + l + 1);
                }
            }
            auto [nodes, weights] = gaussLegendrePoints(mesh(j), mesh(j + 1));

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
        cout << "\n--- N = " << N << " ---" << endl;
        auto start = std::chrono::high_resolution_clock::now();
        ddgic.RKDDG();
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::micro> duration = end - start;
        cout << "Time s = " << duration.count() / 1e6 << " seconds" << endl;
        VectorXd X = VectorXd::LinSpaced(640000 + 1, Xa, Xb);
        h_values.push_back(ddgic.getn(X));
    }
    // --- 误差分析 ---
    cout << "beta0 = " << beta0 << ", beta1 = " << beta1 << ", T = " << T << ", CFL = " << CFL << endl;

    for (size_t i = 0; i < Ns.size() - 1; ++i) {
        
        VectorXd error_vec = h_values[i] - h_values[i + 1];

        double h_i = (Xb - Xa) / Ns[i];
        L2_errors[i] = error_vec.norm(); // 近似L2积分范数
        LInf_errors[i] = error_vec.lpNorm<Infinity>();

        std::cout << "N: " << Ns[i] << " (h=" << h_i << ")"
            << ", L2 Error: " << L2_errors[i]
            << ", LInf Error: " << LInf_errors[i] << std::endl;
    }

    // --- 计算收敛阶 ---
    cout << "\n--- Convergence Order Analysis ---" << endl;
    for (size_t i = 0; i < Ns.size() - 1; ++i) {

        double l2_order = log(L2_errors[i] / L2_errors[i + 1]) / log(2.0);
        double l_inf_order = log(LInf_errors[i] / LInf_errors[i + 1]) / log(2.0);

        std::cout << "Order from N=" << Ns[i] << " to N=" << Ns[i + 1] << ":"
            << " L2 Order = " << l2_order
            << ", L_inf Order = " << l_inf_order << std::endl;
    }

    return 0;
    // 计算误差
}