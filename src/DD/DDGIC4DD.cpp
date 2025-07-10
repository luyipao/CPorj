#include <iostream>
#include <gsl/gsl_integration.h>
#include "../include/gaussLegendre.h"
#include <gsl/gsl_sf_legendre.h>
#include<eigen3/unsupported/Eigen/KroneckerProduct> 
#include "../include/showProgressBar.h"
#include <chrono>
#include "PhyConst.h"
#include "../include/discretize.h"
#include <unordered_map>
#include <Eigen/Sparse>

using namespace Eigen;
using namespace std;

class DDGIC4DD {
private:
    unordered_map<int, pair<VectorXd, VectorXd>> gaussPoints;
    SparseMatrix<double> A_sparse;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> A_ldlt;
    // L() 函数中用到的常量
    const size_t LL = 0, LR = 1, RL = 2, RR = 3;
    VectorXd uxHelp; // 用于存储 ux 的帮助向量
    VectorXd uxxHelp; // 用于存储 uxx 的帮助向量

public:
    int N; // Number of cells
    double Xa; // Left boundary
    double Xb; // Right boundary
    int k; // Polynomial max degree
    double beta0;
    double beta1;
    double T; // Total time
    double CFL; // CFL condition
    double dt; // Time step size
    double h; // Grid size
    MatrixXd C; // Coefficient matrix
    VectorXd mesh; // Mesh points
    MatrixXd n_d_C;
public:
    DDGIC4DD(int N, double Xa, double Xb, int k, double beta0, double beta1, double T, double CFL)
        : N(N), Xa(Xa), Xb(Xb), k(k), beta0(beta0), beta1(beta1), T(T), CFL(CFL) {
        h = (Xb - Xa) / N;
        dt = CFL * h * h; // CFL condition
        mesh = VectorXd::LinSpaced(N + 1, Xa, Xb); // Mesh points
        for (int j = 0; j < N; ++j) {
            gaussPoints[j] = gaussLegendrePoints(mesh(j), mesh(j + 1));
        }

        // 构建质量矩阵Aj
        MatrixXd Aj(k + 1, k + 1);
        for (int l = 0; l <= k; l++) {
            for (int m = 0; m <= k; m++) {
                Aj(l, m) = h / (m + l + 1.0);
            }
        }
        A_sparse.resize((k + 1) * N, (k + 1) * N);
        vector<Eigen::Triplet<double>> tripletList;
        tripletList.reserve(N * (k + 1) * (k + 1));
        for (int j = 0; j < N; ++j) { // 遍历 N 个对角块
            int block_offset = j * (k + 1);
            for (int l = 0; l <= k; ++l) { // 块内行索引
                for (int m = 0; m <= k; ++m) { // 块内列索引
                    double value = h / (m + l + 1.0);
                    int global_row = block_offset + l;
                    int global_col = block_offset + m;
                    tripletList.push_back(Eigen::Triplet<double>(global_row, global_col, value));
                }
            }
        }
        A_sparse.setFromTriplets(tripletList.begin(), tripletList.end());
        A_ldlt.compute(A_sparse);

        uxHelp = VectorXd::LinSpaced(k + 1, 0, k);
        uxxHelp = uxHelp.array() * (uxHelp.array() - 1);
    }
    /**
     * Sets the coefficient matrix of n
     * @param coeff : every column represents the coefficients for a cell.
     * @throws std::invalid_argument if the size of coeff does not match (k + 1) x N.
     */
    void setn(const MatrixXd& coeff) {
        if (coeff.rows() != k + 1 || coeff.cols() != N) {
            throw std::invalid_argument("Coefficient matrix size mismatch.");
        }
        this->C = coeff;
    }
    /**
     * @brief Sets the coefficient matrix of n_d
     *
     * @param coeff : every column represents the coefficients for a cell.
     * @throws std::invalid_argument if the size of coeff does not match (k + 1) x N.
     */
    void setn_d(const MatrixXd& coeff) {
        if (coeff.rows() != k + 1 || coeff.cols() != N) {
            throw std::invalid_argument("Coefficient matrix size mismatch.");
        }
        this->n_d_C = coeff;
    }

    /**
     * @brief 计算 E_x 的原函数（从 Xa 到 X 的定积分）
     *
     * @param X              要求积分的点
     * @return VectorXd      在每个点 X(i) 处的定积分
     */
    VectorXd getPriE_x(const VectorXd& X) {
        const MatrixXd Coeff = C - n_d_C; // 被积函数的多项式系数
        VectorXd Y(X.size());
        Y.setZero();

        // 优化：预计算系数矩阵的累积和，以避免在循环中重复求和。
        // CumSumCoeff(m, j) = sum_{i=0 to j} Coeff(m, i)
        MatrixXd CumSumCoeff(k + 1, N);
        CumSumCoeff.col(0) = Coeff.col(0);
        for (int j = 1; j < N; ++j) {
            CumSumCoeff.col(j) = CumSumCoeff.col(j - 1) + Coeff.col(j);
        }

        // 2. 获取所有点所在的单元格索引 (0-based)
        const VectorXi bins = discretize_uniform(X, mesh);
        for (int i = 0; i < X.size(); ++i) {
            const double xi = X(i);
            if (xi == Xa) {
                Y(i) = 0.0;
                continue;
            }
            const int bin_idx = bins(i);
            if (bin_idx < 0) {
                Y(i) = std::numeric_limits<double>::quiet_NaN();
                continue;
            }
            const double u = (xi - mesh(bin_idx)) / h;
            double integral_at_xi = 0.0;

            // 对每个多项式阶数 m 进行计算
            for (int m = 0; m <= k; ++m) {
                const double m_plus_1 = static_cast<double>(m + 1);
                double full_cells_integral_contrib = 0.0;
                if (bin_idx > 0) {
                    // sum(Coeff(m, 0:bin_idx-1)) * h / (m+1)
                    full_cells_integral_contrib = CumSumCoeff(m, bin_idx - 1) * h / m_plus_1;
                }
                double partial_cell_integral_contrib = Coeff(m, bin_idx) * h * pow(u, m_plus_1) / m_plus_1;
                integral_at_xi += full_cells_integral_contrib + partial_cell_integral_contrib;
            }
            Y(i) = integral_at_xi;
        }
        return Y;
    }

    /**
     * @brief getPriE_x 的标量版本
     */
    double getPriE_x(const double x) {
        VectorXd X_vec(1);
        X_vec(0) = x;
        return getPriE_x(X_vec)(0);
    }

    size_t discretize(const double x) {
        if (x < Xa || x > Xb) {
            throw std::out_of_range("x is out of bounds.");
        }
        if (x == Xb) {
            return N - 1;
        }
        return static_cast<size_t>((x - Xa) / h);
    }

    /**
     * @brief Gets the value of E at a given vector X.
     * Gets the value of E at mesh X.
     * @param x : The point at which to evaluate E.
     * @return The value of E at point x.
     */
    VectorXd getE(const VectorXd& X) {
        auto [nodes, weights] = gaussLegendrePoints(0, 0.6, 10);
        double E0 = 0;
        E0 += getPriE_x(nodes).dot(weights);
        auto [nodes2, weights2] = gaussLegendrePoints(0, 0.4, 10);
        E0 += getPriE_x(nodes2).dot(weights2);
        E0 += 0.4 * getPriE_x(0.6);
        E0 *= 0.001546423010635;
        VectorXd y = -0.001546423010635 * getPriE_x(X).array() + getPriE_x(0) + E0 - 1.5;
        return y;
    }

    /**
     * Gets the value of n at a given point x.
     * @param x : The point at which to evaluate n.
     * @return The value of n at point x.
     */
    double getn(const double x) {
        size_t idx = discretize(x);
        double Cellx = (x - mesh(idx)) / h;
        VectorXd exp(k + 1);
        exp(0) = 1;
        for (int i = 1; i <= k; ++i) {
            exp(i) = exp(i - 1) * Cellx;
        }
        double y = C.col(idx).dot(exp);
        return y;
    }
    /**
     * Gets the value of n at mesh X.
     * @param X : The points at which to evaluate n.
     * @return The values of n at points X.
     */
    VectorXd getn(const VectorXd& X) {
        Index M = X.size();
        VectorXi idx = discretize_uniform(X, mesh);
        // 2. Compute normalized coordinate in each cell for all points (vectorized).
        // If mesh is uniform:
        ArrayXd Xa_array = ArrayXd::Constant(M, Xa);
        ArrayXd leftEdges = Xa_array + idx.cast<double>().array() * h;
        ArrayXd Cellx = (X.array() - leftEdges) / h;

        // 3. Polynomial evaluation using Horner's method (vectorized over points).
        // Iterate over polynomial degree (small loop over k, no loop over data size M).
        VectorXd Y = VectorXd::Zero(M);
        for (int j = k; j >= 0; --j) {
            // Gather the j-th coefficient for each point's cell index into a vector
            VectorXd coeff_j = idx.unaryExpr([&](int cell) { return C(j, cell); });
            // Update Y: Y = Y * Cellx + coeff_j, done as array operations on M elements
            Y.array() = Y.array() * Cellx + coeff_j.array();
        }
        return Y;
    }
    /**
     * Gets the value of n_d at a given point x.
     * @param x : The point at which to evaluate n_d.
     * @return The value of n_d at point x.
     */
    double getn_x(const double x) {
        size_t idx = discretize(x);
        double Cellx = (x - mesh(idx)) / h;
        double y = 0.0;
        double tempx = 1;
        for (int i = 1; i <= k; ++i) {
            y += C(i, idx) * i * tempx;
            tempx *= Cellx;
        }
        return y / h;
    }

    VectorXd getn_x(const VectorXd& X) {
        VectorXd result(X.size());
        for (int i = 0; i < X.size(); ++i) {
            result(i) = getn_x(X(i));
        }
        return result;
    }

    VectorXd phi_x(const VectorXd& X, int k) {
        VectorXd Y(X.size());
        for (int i = 0; i < X.size(); ++i) {
            double xi = X(i);
            size_t idx = discretize(xi);
            double Cellx = (xi - mesh(idx)) / h;
            Y(i) = k / h * pow(Cellx, k - 1);
        }
        return Y;
    }


    MatrixXd L() {

        // 初始化u, ux, uxx
        MatrixXd u = MatrixXd::Zero(4, N);
        MatrixXd ux = MatrixXd::Zero(4, N);
        MatrixXd uxx = MatrixXd::Zero(4, N);

        u.row(LL).tail(N - 1) = C.block(0, 0, k + 1, N - 1).colwise().sum(); // LL
        u.row(RR).head(N - 1) = C.row(0).tail(N - 1); // RR
        u.row(LR) = C.row(0);
        u.row(RL) = C.colwise().sum();

        if (k >= 1) {
            VectorXd C_dot_uxHelp = C.transpose() * uxHelp;;
            ux.row(LL).tail(N - 1) = C_dot_uxHelp.head(N - 1); // LL
            ux.row(RR).head(N - 1) = C.row(1).tail(N - 1); // RR
            ux.row(LR) = C.row(1);
            ux.row(RL) = C_dot_uxHelp;
            ux = ux.array() / h;
        }
        if (k >= 2) {
            VectorXd C_dot_uxxHelp = C.transpose() * uxxHelp;;
            uxx.row(LL).tail(N - 1) = C_dot_uxxHelp.head(N - 1); // LL
            uxx.row(RR).head(N - 1) = 2 * C.row(2).tail(N - 1).array();
            uxx.row(LR) = 2 * C.row(2);
            uxx.row(RL) = C_dot_uxxHelp;
            uxx = uxx.array() / (h * h);
        }

        // 周期边界条件
        if (k >= 0) {
            u(LL, 0) = u(RL, N - 1);
            u(RR, N - 1) = u(LR, 0);
        }
        if (k >= 1) {
            ux(LL, 0) = ux(RL, N - 1);
            ux(RR, N - 1) = ux(LR, 0);
        }
        if (k >= 2) {
            uxx(LL, 0) = uxx(RL, N - 1);
            uxx(RR, N - 1) = uxx(LR, 0);
        }

        // 储存电场E，和高斯积分点上的E，高斯积分点上的phi_x
        /**
        VectorXd E = getE(mesh);
        vector<VectorXd> E_gauss;
        vector<vector<VectorXd>> phi_x_gauss(N, vector<VectorXd>(k + 1));
        VectorXd combineMeshGauss;
        combineMeshGauss << mesh;
        for (int j = 0; j < N; ++j) {
            auto [nodes, weights] = gaussPoints[j];
            combineMeshGauss << nodes;
            for (int i = 0; i <= k; ++i) {
                phi_x_gauss[j][i] = phi_x(nodes, i);
            }
        }
        VectorXd EHelp = getE(combineMeshGauss);

        for (int j = 0; j < N; ++j) {
            E = EHelp.head(mesh.size());
            E_gauss.emplace_back(EHelp.segment(mesh.size() + j * gaussPoints[0].first.size(), gaussPoints[0].first.size()));
        }

*/
// 1. 构建一个包含所有需要计算 E 的点的组合向量
//    预先计算总尺寸，避免动态增长
        const int num_mesh_points = mesh.size();
        const int num_gauss_nodes = gaussPoints[0].first.size(); // 假设所有单元的高斯点数相同
        const int total_points = num_mesh_points + N * num_gauss_nodes;

        VectorXd combined_points(total_points);
        combined_points.head(num_mesh_points) = mesh;
        for (int j = 0; j < N; ++j) {
            const auto& [nodes, weights] = gaussPoints[j];
            combined_points.segment(num_mesh_points + j * num_gauss_nodes, num_gauss_nodes) = nodes;
        }

        VectorXd E_all = getE(combined_points);
        VectorXd E = E_all.head(num_mesh_points);
        vector<VectorXd> E_gauss(N);
        for (int j = 0; j < N; ++j) {
            E_gauss[j] = E_all.segment(num_mesh_points + j * num_gauss_nodes, num_gauss_nodes);
        }

        vector<vector<VectorXd>> phi_x_gauss(N, vector<VectorXd>(k + 1));
        for (int j = 0; j < N; ++j) {
            for (int i = 0; i <= k; ++i) {
                phi_x_gauss[j][i] = phi_x(gaussPoints[j].first, i);
            }
        }
        // 计算通量项B
        MatrixXd B = MatrixXd::Zero(k + 1, N);
        VectorXd tempB1 = PhyConst::mu * (
            E.segment(1, N).cwiseMin(0.0).cwiseProduct(u.row(RL).transpose()) +
            E.segment(1, N).cwiseMax(0.0).cwiseProduct(u.row(RR).transpose())
            );

        VectorXd tempB2 = beta0 * (u.row(RR) - u.row(RL)).transpose() / h;
        VectorXd tempB3 = 0.5 * (ux.row(RR) + ux.row(RL)).transpose();
        VectorXd tempB4 = beta1 * h * (uxx.row(RR) - uxx.row(RL)).transpose();

        for (int l = 0; l <= k; l++) {
            B.row(l) = (tempB1 + PhyConst::tau * PhyConst::theta *
                (tempB2 + tempB3 + tempB4)).transpose();

            if (l == 0) {
                VectorXd tempB5 = PhyConst::mu * (
                    E.segment(0, N).cwiseMin(0.0).cwiseProduct(u.row(LL).transpose()) +
                    E.segment(0, N).cwiseMax(0.0).cwiseProduct(u.row(LR).transpose())
                    );

                VectorXd tempB6 = beta0 * (u.row(LR) - u.row(LL)).transpose() / h;
                VectorXd tempB7 = 0.5 * (ux.row(LR) + ux.row(LL)).transpose();
                VectorXd tempB8 = beta1 * h * (uxx.row(LR) - uxx.row(LL)).transpose();

                B.row(l) -= (tempB5 + PhyConst::tau * PhyConst::theta *
                    (tempB6 + tempB7 + tempB8)).transpose();
            }
        }

        // 计算积分项D
        MatrixXd D = MatrixXd::Zero(k + 1, N);
        for (int j = 0; j < N; ++j) {
            auto [nodes, weights] = gaussPoints[j];
            for (int l = 0; l <= k; ++l) {
                D(l, j) = PhyConst::mu * (weights.array() * E_gauss[j].array() * getn(nodes).array() * phi_x_gauss[j][l].array()).sum();
                D(l, j) += PhyConst::tau * PhyConst::theta * (weights.array() * getn_x(nodes).array() * phi_x_gauss[j][l].array()).sum();
            }
        }

        // 计算界面修正项E
        MatrixXd E_mat = MatrixXd::Zero(k + 1, N);
        if (k >= 1) {
            for (int l = 1; l <= k; l++) {
                for (int j = 0; j < N; j++) {
                    double term1 = (u(RR, j) - u(RL, j)) * l / h * pow(1.0, l - 1);
                    double term2 = (u(LR, j) - u(LL, j)) * l / h * (l == 1);
                    E_mat(l, j) = PhyConst::tau * PhyConst::theta * 0.5 * (term1 + term2);
                }
            }
        }

        // 构建右端向量
        MatrixXd b_mat = B - D - E_mat;
        VectorXd b = Eigen::Map<VectorXd>(b_mat.data(), b_mat.size());

        // 求解线性系统
        VectorXd res_vec = A_ldlt.solve(b);

        // 重塑结果为矩阵
        Eigen::Map<MatrixXd> res(res_vec.data(), k + 1, N);

        return res;
    }

    void RK() {
        MatrixXd coeff = C;
        auto k1 = L();
        k1 = coeff + dt * k1;
        setn(k1);
        auto k2 = L();
        k2 = 3.0 / 4.0 * coeff + 1.0 / 4.0 * k1 + 1.0 / 4.0 * dt * k2;
        setn(k2);
        auto k3 = L();
        k3 = 1.0 / 3.0 * coeff + 2.0 / 3.0 * k2 + 2.0 / 3.0 * dt * k3;
        setn(k3);
    }
    void RKDDG() {
        int count = 0;
        int L = T / dt;
        while (count <= L) {
            count++;
            showProgressBar(L, count);
            RK();
        }
    }
};
