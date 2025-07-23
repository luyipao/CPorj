#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <numeric>

using namespace Eigen;
using namespace std;

// S4(z)过渡函数，z在[0, 1]之间
double S4(double z) {
    if (z <= 0.0) return 0.0;
    if (z >= 1.0) return 1.0;
    // 使用霍纳法则（Horner's method）计算 S4(z) = 70z^9 - 315z^8 + 540z^7 - 420z^6 + 126z^5
    return pow(z, 5) * (126 + z * (-420 + z * (540 + z * (-315 + z * 70))));
}

// S4(z) 的一阶导数: S4'(z) = 630z^4(1-z)^4
double S4_prime(double z) {
    if (z <= 0.0 || z >= 1.0) return 0.0;
    return 630.0 * pow(z, 4) * pow(1 - z, 4);
}

// S4(z) 的二阶导数: S4''(z) = 2520z^3(1-z)^3(1-2z)
double S4_double_prime(double z) {
    if (z <= 0.0 || z >= 1.0) return 0.0;
    return 2520.0 * pow(z, 3) * pow(1 - z, 3) * (1 - 2 * z);
}


// 电子密度分布函数 n_d(x)
// 使用单值输入以简化导数计算
double n_d(double x_val) {
    if (x_val < 0 || x_val > 0.6) {
        throw out_of_range("X value out of bounds: " + to_string(x_val));
    }

    if (x_val >= 0 && x_val <= 0.1) {
        return 5e5;
    } else if (x_val > 0.1 && x_val <= 0.15) {
        double x_norm = (x_val - 0.1) / 0.05;
        double s4_val = S4(x_norm);
        return 5e5 * (1 - s4_val) + 2e3 * s4_val;
    } else if (x_val > 0.15 && x_val <= 0.45) {
        return 2e3;
    } else if (x_val > 0.45 && x_val <= 0.5) {
        double x_norm = (x_val - 0.45) / 0.05;
        double s4_val = S4(x_norm);
        return 2e3 * (1 - s4_val) + 5e5 * s4_val;
    } else { // x_val > 0.5 && x_val <= 0.6
        return 5e5;
    }
}

// n_d 的一阶导数 n_d_x(x)
double n_d_x(double x_val) {
    if (x_val < 0 || x_val > 0.6) {
        throw out_of_range("X value out of bounds for derivative: " + to_string(x_val));
    }
    
    if (x_val > 0.1 && x_val <= 0.15) {
        double h_s = 0.05;
        double x_norm = (x_val - 0.1) / h_s;
        return (2e3 - 5e5) / h_s * S4_prime(x_norm);
    } else if (x_val > 0.45 && x_val <= 0.5) {
        double h_s = 0.05;
        double x_norm = (x_val - 0.45) / h_s;
        return (5e5 - 2e3) / h_s * S4_prime(x_norm);
    } else {
        // 在平坦区域，导数为0
        return 0.0;
    }
}

// n_d 的二阶导数 n_d_xx(x)
double n_d_xx(double x_val) {
    if (x_val < 0 || x_val > 0.6) {
        throw out_of_range("X value out of bounds for 2nd derivative: " + to_string(x_val));
    }
    
    if (x_val > 0.1 && x_val <= 0.15) {
        double h_s = 0.05;
        double x_norm = (x_val - 0.1) / h_s;
        return (2e3 - 5e5) / (h_s * h_s) * S4_double_prime(x_norm);
    } else if (x_val > 0.45 && x_val <= 0.5) {
        double h_s = 0.05;
        double x_norm = (x_val - 0.45) / h_s;
        return (5e5 - 2e3) / (h_s * h_s) * S4_double_prime(x_norm);
    } else {
        // 在平坦区域，二阶导数为0
        return 0.0;
    }
}


int main() {
    // --- 参数定义 ---
    // 系数 beta_1, 可根据需要修改
    const double BETA_1 = 1.0 / 12.0; 
    // 离散化区间的数量 N, 增加 N 会提高精度
    const int N = 80; 
    // 定义域
    const double X_MIN = 0.0;
    const double X_MAX = 0.6;
    // 步长 h
    const double h = (X_MAX - X_MIN) / N;

    // 创建网格点 (N+1个点，形成N个区间)
    VectorXd x_grid = VectorXd::LinSpaced(N + 1, X_MIN, X_MAX);
    
    // --- 计算分子 ---
    // 由于 n_d 是 C4 连续函数, { (n_d)_x } = 0
    // 分子项简化为 (beta_1 / 2 * h * (n_d)_xx)^2
    double numerator_sum = 0.0;
    // 对所有内部节点 j=1, ..., N-1 进行求和
    // 注意: 原公式的求和j=1到N可以解释为对N个区间或N个节点的求和。
    // 这里我们采用对所有N+1个节点求和，这在数值积分中更常见，且对于平滑函数结果差异极小。
    for (int j = 0; j <= N; ++j) {
        double n_d_x_val = n_d_x(x_grid(j));
        double nd_xx_val = n_d_xx(x_grid(j));
        numerator_sum += n_d_x_val + BETA_1 / 2 * nd_xx_val * nd_xx_val;
    }
    double numerator = h * numerator_sum;

    // --- 计算分母 ---
    // 使用梯形法则计算积分: integral |(n_d)_x|^2 dx
    double denominator_integral = 0.0;
    for (int j = 0; j <= N; ++j) {
        double nd_x_val = n_d_x(x_grid(j));
            denominator_integral += nd_x_val * nd_x_val;
    }
    denominator_integral *= h;
    
    // --- 计算最终结果 ---
    if (abs(denominator_integral) < 1e-12) {
        cout << "Error: Denominator is zero or too small." << endl;
        return 1;
    }

    double result = numerator / denominator_integral;

    // --- 输出结果 ---
    cout << "Calculation Parameters:" << endl;
    cout << "  BETA_1 = " << BETA_1 << endl;
    cout << "  N (intervals) = " << N << endl;
    cout << "  h (step size) = " << h << endl;
    cout << "---------------------------------" << endl;
    cout << "Numerator value: " << scientific << numerator << endl;
    cout << "Denominator value: " << scientific << denominator_integral << endl;
    cout << "---------------------------------" << endl;
    cout << "Final Result: " << scientific << result << endl;

    return 0;
}