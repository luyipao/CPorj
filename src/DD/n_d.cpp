#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>

using namespace Eigen;
using namespace std;

// 电子密度分布函数，采用 S4(x) 作为过渡函数
Eigen::VectorXd n_d(const VectorXd& X) {
    VectorXd Y(X.size());
    for (int i = 0; i < X.size(); ++i) {
        double x_val = X[i]; // 使用新变量以避免混淆
        if (x_val < 0 || x_val > 0.6) {
            throw out_of_range("X value out of bounds: " + to_string(x_val));
        }

        if (x_val >= 0 && x_val <= 0.1) {
            Y(i) = 5e5;
        }
        else if (x_val > 0.1 && x_val <= 0.15) {
            double x_norm = (x_val - 0.1) / 0.05; // 归一化到 [0, 1]
            
            // 采用霍纳法则（Horner's method
            double tempy = pow(x_norm, 5) * (126 + x_norm * (-420 + x_norm * (540 + x_norm * (-315 + 70 * x_norm))));
            
            Y(i) = 5e5 * (1 - tempy) + 2e3 * tempy;
        }
        else if (x_val > 0.15 && x_val <= 0.45) {
            Y(i) = 2e3; // 注意：你原来的代码这里写的是 x = 2e3，这应该是 Y(i) = 2e3
        }
        else if (x_val > 0.45 && x_val <= 0.5) {
            double x_norm = (x_val - 0.45) / 0.05; // 归一化到 [0, 1]

            // 同样使用 S_4(x)
            double tempy = pow(x_norm, 5) * (126 + x_norm * (-420 + x_norm * (540 + x_norm * (-315 + 70 * x_norm))));

            Y(i) = 2e3 * (1 - tempy) + 5e5 * tempy;
        }
        else if (x_val > 0.5 && x_val <= 0.6) {
            Y(i) = 5e5;
        }
    }
    return Y;
}
