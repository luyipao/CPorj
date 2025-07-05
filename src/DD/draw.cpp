#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include "n_d.cpp"

using namespace Eigen;


void draw(const MatrixXd& C) {
    // 生成x值：从-10到10的200个点
    const int num_points = 1000;
    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(num_points, 0, 0.6);
    int N = C.cols();
    int k = C.rows() - 1; // 假设C的行数
    double h = 0.6 / N; // 网格大小
    VectorXd mesh = VectorXd::LinSpaced(N + 1, 0, 0.6); // 网格点
    // 计算y值
    Eigen::VectorXd y = VectorXd(num_points);
    y.setZero(); // 初始化y为0
    
    for (int i = 0; i < num_points; ++i) {
        double xi = x(i);
        int idx;
        if (xi == 0.6) {
                idx = N - 1;
        }
        else {
            idx = (xi - 0) / h;
        }
        xi = (xi - mesh(idx)) / h;
        for (int j = 0; j <= k; ++j) {
            y(i) += C(j, idx) * pow(xi, j );
        }
    }
    // 保存数据到文件
    std::ofstream file("data.txt");
    if (file.is_open()) {
        for (int i = 0; i < x.size(); ++i) {
            file << x(i) << " " << y(i) << "\n";
        }
        file.close();
        std::cout << "数据已保存到 data.txt" << std::endl;
    }
    else {
        std::cerr << "无法打开文件！" << std::endl;
    }
}