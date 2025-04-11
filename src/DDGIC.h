#ifndef DDGIC_H
#define DDGIC_H
#include <Eigen/Dense>

using namespace Eigen;

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
};

#endif