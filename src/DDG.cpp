#include <iostream>
#include <gsl/gsl_integration.h>
#include "./include/gaussLegendre.h"
#include <gsl/gsl_sf_legendre.h>
#include<eigen3/unsupported/Eigen/KroneckerProduct> 
#include <chrono>

using namespace std;

double u(double x) {
    return sin(x);
}

/**
 * _Pragma: coeff: u coeff in time t^n
 * k: max k - 1 order polynomial basis, all k base functions
 */
 /**
  * @brief L function for DDG without interface correction
  * @param coeff: u coeff in time t^n
  * @param A: mass matrix
  * @param B: mass matrix
  * @param D: mass matrix
  * @param E: mass matrix
  * @return u_t coeff for problem u_t = L(u)
  */
VectorXd L(const VectorXd& coeff, const MatrixXd& A, const MatrixXd& B, const MatrixXd& D, const MatrixXd& E) {
    MatrixXd M = B - D - E;
    return A.colPivHouseholderQr().solve(M * coeff);
}
/**
 * @brief RK function for DDG without interface correction
 * @param coeff: u coeff in time t^n
 * @param A: mass matrix
 * @param B: mass matrix
 * @param D: mass matrix
 * @param E: mass matrix
 * @param dt: time step
 * @return u coeff in time t^{n+1}
 */
VectorXd RK(const VectorXd& coeff, const MatrixXd& A, const MatrixXd& B, const MatrixXd& D, const MatrixXd& E, const double dt) {
    auto k1 = L(coeff, A, B, D, E);
    k1 = coeff + dt * k1;
    auto k2 = L(k1, A, B, D, E);
    k2 = 3.0 / 4.0 * coeff + 1.0 / 4.0 * k1 + 1.0 / 4.0 * dt * k2;
    auto k3 = L(k2, A, B, D, E);
    k3 = 1.0 / 3.0 * coeff + 2.0 / 3.0 * k2 + 2.0 / 3.0 * dt * k3;
    return k3;
}
/**
 * @brief RKDG function for DDG without interface correction
 * @param Xa: left boundary
 * @param Xb: right boundary
 * @param N: number of cells
 * @param CFL: CFL number
 * @param k: polynomial basis with degree 0 -> k - 1, total k basis functions
 * @param T: Final time
 * @return error vector
 */
VectorXd RKDG(const double Xa, const double Xb, int N, double CFL, int k, double T, double beta0, double beta1) {
    VectorXd mesh(N + 1);
    double meshSize = (Xb - Xa) / N;
    for (int i = 0; i <= N; ++i) {
        mesh(i) = Xa + i * meshSize;
    }

    // construct initial function coeff
    vector<VectorXd> coeff(N, VectorXd::Zero(k)); // 每个单元的系数向量
    for (int j = 0; j < N; ++j) {    // 遍历cell
        // 构建质量矩阵M和右端项b
        MatrixXd M = MatrixXd::Zero(k, k);
        for (int i = 1; i <= k; ++i) {
            for (int j = 1; j <= k; ++j) {
                M(i - 1, j - 1) = meshSize / (i + j - 1);
            }
        }

        VectorXd b = VectorXd::Zero(k);
        for (int i = 1; i <= k; ++i) {
            int gslNum = 10; // 高斯积分点个数
            auto [nodes, weights] = gaussLegendrePoints(mesh(j), mesh(j + 1), gslNum);
            for (int l = 0; l < gslNum; ++l) {
                b(i - 1) += weights(l) * u(nodes(l)) * pow((nodes(l) - mesh(j)) / meshSize, i - 1);
            }
        }

        coeff[j] = M.colPivHouseholderQr().solve(b);
    }

    // construct mass function
    MatrixXd Aj = MatrixXd::Zero(k, k);
    for (int m = 1; m <= k; ++m) {
        for (int l = 1; l <= k; ++l) {
            Aj(l - 1, m - 1) = meshSize / (m + l - 1);
        }
    }
    MatrixXd I = MatrixXd::Identity(N, N);
    MatrixXd A = kroneckerProduct(I, Aj);

    MatrixXd Bj = MatrixXd::Zero(k, k);
    MatrixXd Bjpp = MatrixXd::Zero(k, k);
    for (int l = 1; l <= k; ++l) {
        for (int m = 1; m <= k; ++m) {
            if (m == 1) {
                Bj(l - 1, m - 1) = -beta0 / meshSize;
            }
            else if (m == 2) {
                Bj(l - 1, m - 1) = -beta0 / meshSize + 0.5 * (m - 1) / meshSize;
            }
            else {
                Bj(l - 1, m - 1) = -beta0 / meshSize + 0.5 * (m - 1) / meshSize - beta1 * (m - 1) * (m - 2) / meshSize;
            }
        }
        if (k == 1) {
            Bjpp(l - 1, 0) = beta0 / meshSize;
        }
        else if (k == 2) {
            Bjpp(l - 1, 0) = beta0 / meshSize;
            Bjpp(l - 1, 1) = 0.5 / meshSize;
        }
        else if (k == 3) {
            Bjpp(l - 1, 0) = beta0 / meshSize;
            Bjpp(l - 1, 1) = 0.5 / meshSize;
            Bjpp(l - 1, 2) = beta1 * 2 / meshSize;
        }
    }

    MatrixXd Ip = MatrixXd::Zero(N, N);
    for (int i = 0; i < N; ++i) {
        Ip(i, (i + 1) % N) = 1;  // 周期性右邻居
    }
    MatrixXd B = kroneckerProduct(I, Bj) + kroneckerProduct(Ip, Bjpp);

    MatrixXd Dj = MatrixXd::Zero(k, k);
    MatrixXd Djss = MatrixXd::Zero(k, k);
    for (int l = 1; l == 1; ++l) {
        for (int m = 1; m <= k; ++m) {
            if (m == 1) {
                Djss(l - 1, m - 1) = -beta0 / meshSize;
            }
            else if (m == 2) {
                Djss(l - 1, m - 1) = -beta0 / meshSize + 0.5 * (m - 1) / meshSize;
            }
            else {
                Djss(l - 1, m - 1) = -beta0 / meshSize + 0.5 * (m - 1) / meshSize - beta1 * (m - 1) * (m - 2) / meshSize;
            }
        }
        if (k == 1) {
            Dj(l - 1, 0) = beta0 / meshSize;
        }
        else if (k == 2) {
            Dj(l - 1, 0) = beta0 / meshSize;
            Dj(l - 1, 1) = 0.5 / meshSize;
        }
        else if (k == 3) {
            Dj(l - 1, 0) = beta0 / meshSize;
            Dj(l - 1, 1) = 0.5 / meshSize;
            Dj(l - 1, 2) = beta1 * 2 / meshSize;
        }
    }

    MatrixXd DU = MatrixXd::Zero(N, N);
    for (int i = 0; i < N; ++i) {
        DU(i, (i - 1 + N) % N) = 1; // 周期性左邻居
    }
    MatrixXd D = kroneckerProduct(I, Dj) + kroneckerProduct(DU, Djss);

    MatrixXd Ej = MatrixXd::Zero(k, k);
    for (int l = 2; l <= k; ++l) {
        for (int m = 2; m <= k; ++m) {
            Ej(l - 1, m - 1) = (m - 1) * (l - 1) / meshSize / (m + l - 3);
        }
    }
    MatrixXd E = kroneckerProduct(I, Ej);

    VectorXd C = VectorXd::Zero(N * k);
    for (int j = 0; j < N; ++j) {
        for (int m = 1; m <= k; ++m) {
            C(j * k + m - 1) = coeff[j](m - 1);
        }
    }

    // RK 迭代
    double dt = CFL * meshSize * meshSize;
    int count = 0;
    while (count * dt <= T) {
        count++;
        C = RK(C, A, B, D, E, dt);
    }
    // 将结果重新导入系数
    for (int j = 0; j < N; ++j) {
        for (int m = 1; m <= k; ++m) {
            coeff[j](m - 1) = C(j * k + m - 1);
        }
    }
    // 计算误差
    int sampleNum = 1000;
    VectorXd error(sampleNum);
    for (int j = 0; j < sampleNum; ++j) {
        double x = (Xb - Xa) / (sampleNum)*j;
        double result = 0;
        auto it = lower_bound(mesh.data(), mesh.data() + mesh.size(), x);
        int p = max(0, int(it - mesh.data()) - 1);
        for (int i = 1; i <= k; ++i) {
            result += coeff[p](i - 1) * pow((x - mesh(p)) / meshSize, i - 1);
        }
        error(j) = abs(result - exp(-T) * sin(x));
    }
    return error;
}

int main() {
    auto start = std::chrono::system_clock::now();
    double CFL = 0.0005;
    double T = 1;
    double Xa = 0;
    double Xb = 2 * M_PI;
    int k = 4; // k is the polynomial degree, k-1 is the max order of polynomial basis, i.e., you should plus 1 when using k 
    double theta = 0.5;
    double beta0 = 7 / (12 * theta);
    double beta1 = theta / 6;
    vector<int> Nlist = { 10 };
    vector<double> errorInf;
    vector<double> errorL2;
    for (auto N : Nlist) {
        VectorXd error = RKDG(Xa, Xb, N, CFL, k, T, beta0, beta1);
        double eInf = 0;
        double eL2 = 0;
        for (int i = 0; i < error.size(); ++i) {
            eInf = max(eInf, abs(error(i)));
            eL2 += error(i) * error(i);
        }
        errorInf.emplace_back(eInf);
        errorL2.emplace_back(sqrt(eL2 / (1000)));
    }
    cout << "误差分析:" << endl;
    cout << "N\tL2 error\tInf error" << endl;
    for (int i = 0; i < errorInf.size(); ++i) {
        cout << "N=" << Nlist[i] << ", ErrorInf: " << errorInf[i] << endl;
    }
    for (int i = 1; i < errorInf.size(); ++i) {
        double rate = log(errorInf[i - 1] / errorInf[i]) / log(2.0);
        cout << "收敛率 (N=" << Nlist[i - 1] << "→" << Nlist[i] << "): " << rate << endl;
    }
    for (int i = 0; i < errorL2.size(); ++i) {
        cout << "N=" << Nlist[i] << ", L2 errors: " << errorL2[i] << endl;
    }
    for (int i = 1; i < errorL2.size(); ++i) {
        double rate = log(errorL2[i - 1] / errorL2[i]) / log(2.0);
        cout << "收敛率 (N=" << Nlist[i - 1] << "→" << Nlist[i] << "): " << rate << endl;
    }
    auto end = std::chrono::system_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
    return 0;
}

