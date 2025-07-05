#include <iostream>
#include <gsl/gsl_integration.h>
#include "./include/gaussLegendre.h"
#include <gsl/gsl_sf_legendre.h>
#include<eigen3/unsupported/Eigen/KroneckerProduct> 
#include "./include/showProgressBar.h"
#include <chrono>
#include "DDGIC.h"

using namespace Eigen;

double u(double x) {
    return exp(x);
}
double lbdy(double x) {
    return exp(x);
}
double rbdy(double x) {
    return exp(x);
}
double lbdyDiff(double x) {
    return exp(x);
}
double rbdyDiff(double x) {
    return exp(x);
}
double lbdyDD(double x) {
    return exp(x);
}
double rbdyDD(double x) {
    return exp(x);
}
// 目前来看，第一个L运算后就得到发散的结果可能是A和b错了，也可能是A条件数太大了
MatrixXd L(DDGIC& ddgic) {
    int LL = 0, LR = 1, RL = 2, RR = 3;
    int k = ddgic.k;
    int N = ddgic.N;
    double beta0 = ddgic.beta0;
    double beta1 = ddgic.beta1;
    double h = ddgic.h;
    MatrixXd C = ddgic.C;

    // 计算质量矩阵
    MatrixXd Aj(k + 1, k + 1);
    for (int l = 0; l <= k; ++l) {
        for (int m = 0; m <= k; ++m) {
            Aj(l, m) = h / (l + m + 1);
        }
    }
    MatrixXd I4A = MatrixXd::Identity(N, N);
    MatrixXd A = kroneckerProduct(I4A, Aj);

    JacobiSVD<MatrixXd> svd(A);
    double cond = svd.singularValues()(0)
        / svd.singularValues()(svd.singularValues().size() - 1);
    // 定义解函数及其导数在边界的左右极限

    MatrixXd u(N, 4);
    MatrixXd ux = MatrixXd::Zero(N, 4);
    MatrixXd uxx = MatrixXd::Zero(N, 4);
    VectorXd natureNums = VectorXd::LinSpaced(k + 1, 0, k);
    VectorXd uxHelp = natureNums;
    VectorXd uxxHelp = natureNums.array() * (natureNums.array() - 1);

    for (int j = 0; j < N; ++j) {
        if (j == 0) {
            u(0, LL) = lbdy(0);
        }
        else {
            u(j, LL) = C.col(j - 1).sum();
        }
        u(j, 1) = C.col(j)(0);
        u(j, 2) = C.col(j).sum();
        if (j == N - 1) {
            u(j, 3) = rbdy(1);
        }
        else {
            u(j, 3) = C.col(j + 1)(0);
        }
    }

    for (int j = 0; j < N; ++j) {
        if (j != 0) {
            ux(j, LL) = C.col(j - 1).dot(uxHelp);
        }
        if (k >= 1) {
            ux(j, LR) = C.col(j)(1);
        }
        ux(j, RL) = C.col(j).dot(uxHelp);
        if (j != N - 1 && k >= 1) {
            ux(j, RR) = C.col(j + 1)(1);
        }
    }
    ux /= h;
    ux(0, 0) = lbdyDiff(0);
    ux(N - 1, 3) = rbdyDiff(1);

    for (int j = 0; j < N; ++j) {
        if (j != 0) {
            uxx(j, LL) = C.col(j - 1).dot(uxxHelp);
        }
        if (k >= 2) {
            uxx(j, LR) = 2 * C.col(j)(2);
        }
        uxx(j, RL) = C.col(j).dot(uxxHelp);
        if (j != N - 1 && k >= 2) {
            uxx(j, RR) = 2 * C.col(j + 1)(2);
        }
    }
    uxx /= h * h;
    uxx(0, 0) = lbdyDD(0);  // 处理dirichlet边界条件
    uxx(N - 1, 3) = rbdyDD(1);

    MatrixXd uu = u.array() * u.array();
    MatrixXd uux = u.array() * ux.array();
    MatrixXd uxux = ux.array() * ux.array();
    MatrixXd uuxx = u.array() * uxx.array();

    MatrixXd B = MatrixXd::Zero(k + 1, N);
    for (int l = 0; l <= k; ++l) {
        VectorXd BlR = VectorXd::Zero(k + 1);
        VectorXd BlL = VectorXd::Zero(k + 1);
        VectorXd theta = VectorXd::Zero(k + 1);

        theta = u.col(RR).cwiseAbs().cwiseMax(u.col(RL).cwiseAbs());
        VectorXd temp1 = (theta.array() * (u.col(RR) - u.col(RL)).array());
        BlR = 0.5 * uu.col(RR) + 0.5 * uu.col(RL);
        BlR -= temp1;
        theta = u.col(LR).cwiseAbs().cwiseMax(u.col(LL).cwiseAbs());
        BlL = 0.5 * uu.col(LR) + 0.5 * uu.col(LL) - (VectorXd)(theta.array() * (u.col(LR) - u.col(LL)).array());

        if (l == 0) {
            B.row(l) = BlR - BlL;
        }
        else {
            B.row(l) = BlR;
        }
    }
    B = B.reshaped() / 2.0;


    MatrixXd D = MatrixXd::Zero(k + 1, N);
    for (int l = 0; l <= k; ++l) {
        VectorXd DlR = VectorXd::Zero(N);
        VectorXd DlL = VectorXd::Zero(N);
        VectorXd temp1, temp2, temp3;
        if (l == 0) {
            DlR = beta0 / h * (uu.col(RR) - uu.col(RL)) + (uux.col(RR) + uux.col(RL)) + beta1 * h * 2 * (uxux.col(RR) - uxux.col(RL) + uuxx.col(RR) - uuxx.col(RL));
            DlL = beta0 / h * (uu.col(LR) - uu.col(LL)) + (uux.col(LR) + uux.col(LL)) + beta1 * h * 2 * (uxux.col(LR) - uxux.col(LL) + uuxx.col(LR) - uuxx.col(LL));
        }
        else {
            DlR = beta0 / h * (uu.col(RR) - uu.col(RL)) + (uux.col(RR) + uux.col(RL)) + beta1 * h * 2 * (uxux.col(RR) - uxux.col(RL) + uuxx.col(RR) - uuxx.col(RL));
        }
        D.row(l) = DlR - DlL;
    }
    D = D.reshaped() / 4.0;


    MatrixXd E = MatrixXd::Zero(k + 1, N);
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
    E = E.reshaped() / 2.0;

    MatrixXd F = MatrixXd::Zero(k + 1, N);
    for (int l = 0; l <= k; ++l) {
        if (l == 0) {
            ;
        }
        else {
            F.row(l) = (uu.col(RR) - uu.col(RL)) * l / h + (uu.col(LR) - uu.col(LL)) * l * (l == 1) / h;
        }
    }
    F = F.reshaped() / 8.0;

    VectorXd b = -B + D + E - F;
    return A.colPivHouseholderQr().solve(b).reshaped(k + 1, N);
}

void RK(DDGIC& ddgic) {
    int k = ddgic.k;
    int N = ddgic.N;
    double dt = ddgic.dt;

    MatrixXd coeff = ddgic.C;
    auto k1 = L(ddgic);
    k1 = coeff + dt * k1;
    ddgic.setCoeff(k1);
    auto k2 = L(ddgic);
    k2 = 3.0 / 4.0 * coeff + 1.0 / 4.0 * k1 + 1.0 / 4.0 * dt * k2;
    ddgic.setCoeff(k2);
    auto k3 = L(ddgic);
    k3 = 1.0 / 3.0 * coeff + 2.0 / 3.0 * k2 + 2.0 / 3.0 * dt * k3;
    ddgic.setCoeff(k3);
}
void RKDG(DDGIC& ddgic) {
    double dt = ddgic.dt;
    int k = ddgic.k;
    int N = ddgic.N;
    double T = ddgic.T;
    int count = 0;
    while (count * dt <= T) {
        count++;
        showProgressBar(T / dt, count);
        RK(ddgic);
    }
}

int main() {
    int k = 0;
    int N = 10;
    double T = 0.5;
    double Xa = 0;
    double Xb = 1;
    double beta0 = 2.0;
    double beta1 = 1.0 / 12.0;
    double CFL = 0.0005;
    DDGIC ddgic(N, Xa, Xb, k, beta0, beta1, T, CFL);
    int LL = 0, LR = 1, RL = 2, RR = 3;
    double h = (Xb - Xa) / N; // 网格大小

    VectorXd mesh = VectorXd::LinSpaced(N + 1, Xa, Xb); // mesh points

    // construct initial function coeff
    vector<VectorXd> coeff(N, VectorXd::Zero(k + 1)); // 每个单元的系数向量
    for (int j = 0; j < N; ++j) {    // 遍历cell
        // 构建质量矩阵M和右端项b
        MatrixXd M = MatrixXd::Zero(k + 1, k + 1);
        for (int i = 0; i <= k; ++i) {
            for (int j = 0; j <= k; ++j) {
                M(i, j) = h / (i + j + 1);
            }
        }

        VectorXd b = VectorXd::Zero(k + 1);
        for (int i = 0; i <= k; ++i) {
            int gslNum = 10; // 高斯积分点个数
            auto [nodes, weights] = gaussLegendrePoints(mesh(j), mesh(j + 1), gslNum);
            for (int l = 0; l < gslNum; ++l) {
                b(i) += weights(l) * u(nodes(l)) * pow((nodes(l) - mesh(j)) / h, i);
            }
        }
        coeff[j] = M.colPivHouseholderQr().solve(b);
    }

    MatrixXd C(k + 1, N); // 每列代表一个cell上的系数向量
    for (int j = 0; j < N; ++j) {
        C.col(j) = coeff[j]; // Assign each coeff[j] (which is a VectorXd) as a column of C
    }





    ddgic.setCoeff(C);
   // RKDG(ddgic);
   ddgic.RKDDG();





    for (int j = 0; j < N; ++j) {
        coeff[j] = ddgic.C.col(j);
    }
    // 计算误差
    int sampleNum = 1000;
    VectorXd error(sampleNum);
    for (int j = 0; j < sampleNum; ++j) {
        double x = Xa + j * (Xb - Xa) / (sampleNum - 1);
        double result = 0;
        auto it = lower_bound(mesh.data(), mesh.data() + mesh.size(), x);
        int p = max(0, int(it - mesh.data()) - 1);
        for (int i = 0; i <= k; ++i) {
            result += coeff[p](i) * pow((x - mesh(p)) / h, i);
        }
        error(j) = result - exp(x);
    }
    double eInf = 0;
    for (int i = 0; i < error.size(); ++i) {
        eInf = max(eInf, abs(error(i)));
    }

    cout << (double)eInf << endl;
    return 0;
}


