#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>

using namespace Eigen;
using namespace std;

// 电子密度分布函数 (固定 n=3)
Eigen::VectorXd n_d(const VectorXd& X) {
    VectorXd Y(X.size());
    for (int i = 0; i < X.size(); ++i) {
        double x = X[i];
        if (x < 0 || x > 0.6) {
            throw out_of_range("X value out of bounds: " + to_string(x));
        }
        if (x >= 0 && x <= 0.1) {
            Y(i) = 5e5;
        }
        else if (x > 0.1 && x <= 0.15) {
            x = (x - 0.1) / 0.05;
            double tempy = x * x * x * (6 * x * x - 15 * x + 10);
            Y(i) = 5e5 * (1 - tempy) + 2e3 * tempy;
        }
        else if (x > 0.15 && x <= 0.45) {
            x = 2e3;
        }
        else if (x > 0.45 && x <= 0.5) {
            x = (x - 0.45) / 0.05;
            double tempy = x * x * x * (6 * x * x - 15 * x + 10);
            Y(i) = 2e3 * (1 - tempy) + 5e5 * tempy;
        }
        else if (x > 0.5 && x <= 0.6) {
            Y(i) = 5e5;
        }
    }
    return Y;
}