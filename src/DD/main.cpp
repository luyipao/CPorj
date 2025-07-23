#include <iostream>
#include "draw.cpp"
#include "../include/gaussLegendre.h"
#include "./DDGIC4DD.cpp"
#include <chrono>       // 用于获取时间
#include <ctime>        // 用于格式化时间
#include <iomanip>      // 用于格式化输出 (例如 setprecision)
using namespace Eigen;
using namespace std;

int main() {
    auto start = std::chrono::system_clock::now();
    std::time_t start_time = std::chrono::system_clock::to_time_t(start);
    int k = 2;
    double Xa = 0.00;
    double Xb = 0.60;

    // double beta0 = (k + 1) * (k + 1);
    // double beta1 = 1 / (k + 1.0) / (2 * k);
    double beta0 = 32;
    double beta1 = 1.0 / 12.0;
    double T = 1;
    double CFL;
    if (k == 0) {
        CFL = 0.1;
    }
    else if (k == 1) {
        CFL = 0.1;
    }
    else {
        CFL = 0.1;
    }
    // construct initial function coeff
    vector<int> Ns = {80, 160, 320}; // 网格数
    vector<vector<VectorXd>> Sols;
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
        // ddgic.dt = 0.0000225;
        ddgic.setn(n_C);
        ddgic.setn_d(n_C);
        // cout << "--- N = " << N << " ---" << endl;
        // auto start = std::chrono::high_resolution_clock::now();
        ddgic.RKDDG();
        // auto end = std::chrono::high_resolution_clock::now();
        // std::chrono::duration<double, std::micro> duration = end - start;
        // cout << "Time s = " << duration.count() / 1e6 << " seconds" << endl;
        // cout << "T = " << ddgic.T << ", dt = " << ddgic.dt << endl;
        VectorXd X = VectorXd::LinSpaced(10000 + 1, Xa, Xb);
        VectorXd Y = ddgic.getn(X);
        h_values.push_back(Y);
    }
    // --- 误差分析 ---
    const std::string filename = "ErrorAnalysis.txt";
    std::ofstream result_file(filename, std::ios_base::app);
   if (!result_file.is_open()) {
        std::cerr << "错误: 无法打开结果文件 " << filename << std::endl;
        return 1; // 返回错误码
    }
    

    auto end = std::chrono::system_clock::now();
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);
    // 使用 std::put_time (C++11) 或 strftime 来格式化时间
    result_file << "======================================================================\n";
    result_file << " 开始于: " << std::ctime(&start_time) << " 结束于：" << std::ctime(&end_time); // ctime 会自带换行符
    result_file << "======================================================================\n";

    result_file << std::fixed << std::setprecision(6);

    result_file << "k = " << k << ", beta0 = " << beta0 << ", beta1 = " << beta1 
                << ", T = " << T << ", CFL = " << CFL << std::endl << std::endl;

    for (size_t i = 0; i < Ns.size() - 1; ++i) {
        Eigen::VectorXd error_vec = h_values[i] - h_values[i + 1];
        double h_i = (Xb - Xa) / Ns[i];
        L2_errors[i] = error_vec.norm() * sqrt((Xb - Xa) / error_vec.size()); 
        LInf_errors[i] = error_vec.lpNorm<Eigen::Infinity>();

        result_file << "N: " << Ns[i] << " (h=" << h_i << ")"
                    << ", L2 Error: "  << L2_errors[i]
                    << ", LInf Error: " << LInf_errors[i] << std::endl;
    }

    // --- 计算收敛阶 ---
    result_file << "\n--- 收敛阶 ---\n";
    for (int i = 0; i < static_cast<int>(Ns.size()) - 2; ++i) { // 注意循环边界，避免访问越界
        double l2_order = log(L2_errors[i] / L2_errors[i + 1]) / log(2.0);
        double l_inf_order = log(LInf_errors[i] / LInf_errors[i + 1]) / log(2.0);

        result_file << "Order from N=" << Ns[i] << " to N=" << Ns[i+1] << ":"
                    << " L2 Order = " << l2_order
                    << ", L_inf Order = " << l_inf_order << std::endl;
    }
    result_file << "\n\n"; 
    result_file.close();
    std::cout << "计算完成，结果已追加到 " << filename << std::endl;
    std::chrono::duration<double> elapsed_seconds = end - start;

    // 5. 输出结果
    std::cout << "运行时间: " << elapsed_seconds.count() << "s" << std::endl;

    return 0;
}