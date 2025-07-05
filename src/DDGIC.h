#ifndef DDGIC_H
#define DDGIC_H
#include <iostream>
#include <gsl/gsl_integration.h>
#include "./include/gaussLegendre.h"
#include <gsl/gsl_sf_legendre.h>
#include<eigen3/unsupported/Eigen/KroneckerProduct> 
#include "./include/showProgressBar.h"
#include <chrono>
#include "DDGIC.h"

using namespace Eigen;
using namespace std;
class DDGIC {
public:
    int N;
    double Xa;
    double Xb;
    int k;
    double beta0;
    double beta1;
    double T;
    double CFL;
    double dt;
    double h; // 网格大小
    MatrixXd C; // 系数矩阵


    DDGIC() : N(0), Xa(0.0), Xb(0.0), k(0), beta0(0.0), beta1(0.0), T(0.0), CFL(0.0), h(0.0), dt(0.0) {}

    DDGIC(int N, double Xa, double Xb, int k, double beta0, double beta1, double T, double CFL)
        : N(N), Xa(Xa), Xb(Xb), k(k), beta0(beta0), beta1(beta1), T(T), CFL(CFL) {
        h = (Xb - Xa) / N;
        dt = CFL * h * h; // CFL condition
    }

    void setCoeff(const MatrixXd& coeff) {
        if (coeff.rows() != k + 1 || coeff.cols() != N) {
            throw std::invalid_argument("Coefficient matrix size mismatch.");
        }
        this->C = coeff;
    }

    MatrixXd L();
    void RK();
    void RKDDG();
};

// L 绝对错误
MatrixXd DDGIC::L() {
    // 初始参数定义
    int LL = 0, LR = 1, RL = 2, RR = 3;

    // 计算u的值在单元边界的左右极限
    MatrixXd u(4, N);
    MatrixXd ux(4, N);
    MatrixXd uxx(4, N);

    // Dirichlet 边界条件
    u(LL, 0) = 1;
    ux(LL, 0) = 1;
    uxx(LL, 0) = 1;
    u(RR, N - 1) = exp(1);
    ux(RR, N - 1) = exp(1);
    uxx(RR, N - 1) = exp(1);

    VectorXd uxHelp = VectorXd::LinSpaced(k + 1, 0, k);
    VectorXd uxxHelp = uxHelp.array() * (uxHelp.array() - 1);
    for (int j = 0; j < N; ++j) {
        if (j == 0) {
            u(LR, j) = C(0, j);
            u(RL, j) = C.col(j).sum();
            u(RR, j) = C(0, j + 1);

            if (k >= 1) {
                ux(LR, j) = C(1, j) / h;
                ux(RL, j) = C.col(j).dot(uxHelp) / h;
                ux(RR, j) = C(1, j + 1) / h;
            }

            if (k >= 2) {
                uxx(LR, j) = 2 / (h * h) * C(2, j);
                uxx(RL, j) = 1 / (h * h) * C.col(j).dot(uxxHelp);
                uxx(RR, j) = 2 / (h * h) * C(2, j);
            }
        }
        else if (j == N - 1) {
            u(LL, j) = C.col(j - 1).sum();
            u(LR, j) = C(0, j);
            u(RL, j) = C.col(j).sum();

            if (k >= 1) {
                ux(LL, j) = 1 / h * C.col(j - 1).dot(uxHelp);
                ux(LR, j) = 1 / h * C(1, j);
                ux(RL, j) = 1 / h * C.col(j).dot(uxHelp);
            }

            if (k >= 2) {
                uxx(LL, j) = 1 / (h * h) * C.col(j - 1).dot(uxxHelp);
                uxx(LR, j) = 2 / (h * h) * C(2, j);
                uxx(RL, j) = 1 / (h * h) * C.col(j).dot(uxxHelp);
            }
        }
        else {
            u(LL, j) = C.col(j - 1).sum();
            u(LR, j) = C(0, j);
            u(RL, j) = C.col(j).sum();
            u(RR, j) = C(0, j + 1);

            if (k >= 1) {
                ux(LL, j) = 1 / h * C.col(j - 1).dot(uxHelp);
                ux(LR, j) = 1 / h * C(1, j);
                ux(RL, j) = 1 / h * C.col(j).dot(uxHelp);
                ux(RR, j) = 1 / h * C(1, j + 1);
            }

            if (k >= 2) {
                uxx(LL, j) = 1 / (h * h) * C.col(j - 1).dot(uxxHelp);
                uxx(LR, j) = 2 / (h * h) * C(2, j);
                uxx(RL, j) = 1 / (h * h) * C.col(j).dot(uxxHelp);
                uxx(RR, j) = 2 / (h * h) * C(2, j + 1);
            }
        }
    }
    // 构造u的一些二级矩阵
    MatrixXd uu = u.array() * u.array();
    MatrixXd uux = u.array() * ux.array();
    MatrixXd uxux = ux.array() * ux.array();
    MatrixXd uuxx = u.array() * uxx.array();

    // 构造质量矩阵
    MatrixXd Aj(k + 1, k + 1);
    for (int l = 0; l <= k; ++l) {
        for (int m = 0; m <= k; ++m) {
            Aj(l, m) = h / (m + l + 1);
        }
    }
    MatrixXd IN = MatrixXd::Identity(N, N);
    MatrixXd A = kroneckerProduct(IN, Aj);

    MatrixXd B(k + 1, N);
    RowVectorXd theta(N), tempB(N);
    for (int l = 0; l <= k; ++l) {
        // 打印维度信息
        theta = u.row(RR).cwiseAbs().cwiseMax(u.row(RL).cwiseAbs());
        tempB = theta.array() * (u.row(RR).array() - u.row(RL).array());
        B.row(l) = 0.5 * (uu.row(RR) + uu.row(RL)) - tempB;

        if (l == 0) {
            theta = u.row(LR).cwiseAbs().cwiseMax(u.row(LL).cwiseAbs());
            tempB = theta.array() * (u.row(LR).array() - u.row(LL).array());
            B.row(l) -= 0.5 * (uu.row(LR) + uu.row(LL)) - tempB;
        }
    }
    B /= 2;

    MatrixXd D(k + 1, N);
    RowVectorXd tempD1(N), tempD2(N), tempD3(N);
    for (int l = 0; l <= k; ++l) {
        tempD1 = beta0 / h * (uu.row(RR) - uu.row(RL));
        tempD2 = 0.5 * (uux.row(RR) + uux.row(RL));
        tempD3 = beta1 * h * 2 * (uuxx.row(RR) - uuxx.row(RL) + uxux.row(RR) - uxux.row(RL));
        D.row(l) = tempD1 + tempD2 + tempD3;
        if (l == 0) {
            tempD1 = beta0 / h * (uu.row(LR) - uu.row(LL));
            tempD2 = 0.5 * (uux.row(LR) + uux.row(LL));
            tempD3 = beta1 * h * 2 * (uuxx.row(LR) - uuxx.row(LL) + uxux.row(LR) - uxux.row(LL));
            D.row(l) -= tempD1 + tempD2 + tempD3;
        }
    }
    D /= 4;

    MatrixXd E(k + 1, N);
    vector<MatrixXd> M1(k + 1, MatrixXd::Zero((k + 1), (k + 1)));
    vector<MatrixXd> M2(k + 1, MatrixXd::Zero((k + 1), (k + 1)));
    for (int l = 0; l <= k; ++l) {
        for (int m = 0; m <= k; ++m) {
            for (int n = 0; n <= k; ++n) {
                M1[l](m, n) = (l > 0) ? l * 1.0 / (m + n + l) : 0;
                M2[l](m, n) = (n > 0 && l > 0) ? n * l / (m + n + l - 1) / h : 0;
            }
        }
    }
    for (int j = 0; j < N; ++j) {
        RowVectorXd Cjrow = C.col(j);
        for (int l = 0; l <= k; ++l) {
            double term1 = Cjrow * M1[l] * C.col(j);  // 1x1 矩阵转标量
            double term2 = Cjrow * M2[l] * C.col(j);
            E(l, j) = term1 - term2;
        }
    }
    E /= 2;

    MatrixXd F(k + 1, N);
    for (int l = 1; l <= k; ++l) {
        F.row(l) = (uu.row(RR) - uu.row(RL)) * l / h + (uu.row(LR) - uu.row(LL)) * l * (l == 1) / h;
    }
    F /= 8;

    // Ax + B - D - E + F = 0
    VectorXd b = (-B + D + E - F).reshaped();
    MatrixXd res = A.fullPivHouseholderQr().solve(b);
    res = res.reshaped(k + 1, N);
    return res;
}

void DDGIC::RK() {
    MatrixXd Clast = C;
    MatrixXd k1 = L();
    MatrixXd y1 = Clast + dt * k1;
    setCoeff(y1);
    MatrixXd k2 = L();
    MatrixXd y2 = 3.0 / 4.0 * Clast + 1.0 / 4.0 * k1 + 1.0 / 4.0 * dt * k2;
    setCoeff(y2);
    MatrixXd k3 = L();
    MatrixXd y3 = 1.0 / 3.0 * Clast + 2.0 / 3.0 * k2 + 2.0 / 3.0 * dt * k3;
    setCoeff(y3);
}

void DDGIC::RKDDG() {
    int count = 0;
    while (dt * count <= T) {
        showProgressBar(T / dt, count);
        RK();
        count++;
    }
}

#endif