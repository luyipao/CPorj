
#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>
#include <fstream>
#include <iostream>

using namespace Eigen;
using namespace std;

double smooth_step(double x) {
    return 70 * pow(x, 9) - 315 * pow(x, 8) + 540 * pow(x, 7) - 420 * pow(x, 6) + 126 * pow(x, 5);
}

VectorXd smooth_step(const VectorXd& X) {
    VectorXd Y = VectorXd::Zero(X.size());
    Y = 70 * X.array().pow(9) - 315 * X.array().pow(8) + 540 * X.array().pow(7) - 420 * X.array().pow(6) + 126 * X.array().pow(5);
    return Y;
}

VectorXd n_d(const VectorXd& X) {
    VectorXd Y = VectorXd::Zero(X.size());
    for (int i = 0; i < X.size(); ++i) {
        double x = X[i];
        if (x < 0 || x > 0.6) {
            throw out_of_range("X value out of bounds: ");
        }
        else if (x >= 0 && x <= 0.1) {
            Y(i) = 5e5;
        }
        else if (x > 0.1 && x < 0.15) {
            double x_norm = (x - 0.1) / 0.05;
            Y(i) = 5e5 * (1 - smooth_step(x_norm)) + 2e3 * smooth_step(x_norm);
        }
        else if (x >= 0.15 && x <= 0.45) {
            Y(i) = 2e3;
        }
        else if (x > 0.45 && x < 0.5) {
            double x_norm = (x - 0.45) / 0.05;
            Y(i) = 2e3 * (1 - smooth_step(x_norm)) + 5e5 * smooth_step(x_norm);
        }
        else if (x >= 0.5 && x <= 0.6) {
            Y(i) = 5e5;
        }
    }
    return Y;
}

// g++ n_d.cpp ./../functions/*.cpp -O3 -o test  -lgslcblas -lm -lgsl
// int main() {
//     auto X = VectorXd::LinSpaced(10000, 0, 0.6);
//     VectorXd Y = n_d(X);

//     std::ofstream file("data.txt");
//     if (file.is_open()) {
//         for (int i = 0; i < X.size(); ++i) {
//             file << X(i) << " " << Y(i) << "\n";
//         }
//         file.close();
//         std::cout << "数据已保存到 data.txt" << std::endl;
//     }
//     else {
//         std::cerr << "无法打开文件！" << std::endl;
//     }
// }