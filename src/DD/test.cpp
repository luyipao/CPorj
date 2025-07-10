#include <iostream>
#include <Eigen/Dense>
#include "../include/discretize.h"
#include "n_d.cpp"
#include <fstream>
using namespace std;
using namespace Eigen;

int main() {
    // 创建一个0，1之间随机向量
    VectorXd X = VectorXd::LinSpaced(1000, 0, 0.6); // 创建等距边界
    VectorXd Y = n_d(X);
        // 保存数据到文件
    std::ofstream file("data.txt");
    if (file.is_open()) {
        for (int i = 0; i < X.size(); ++i) {
            file << X(i) << " " << Y(i) << "\n";
        }
        file.close();
        std::cout << "数据已保存到 data.txt" << std::endl;
    }
    else {
        std::cerr << "无法打开文件！" << std::endl;
    }
}