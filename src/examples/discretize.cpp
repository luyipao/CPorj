#include <iostream>
#include <Eigen/Dense>
#include "../include/discretize.h"
using namespace std;
using namespace Eigen;

int main() {
    // 创建一个0，1之间随机向量
    VectorXd vec = VectorXd::Random(10);
    vec = vec.array() * 0.5 + 0.5;
    VectorXd edges = VectorXd::LinSpaced(11, 0, 1); // 创建等距边界
    vec[0] = 0;
    vec[1] = 1;
    VectorXi indices = discretize_uniform(vec, edges);

    for (int i = 0; i < indices.size(); ++i) {
        cout << "vec[" << i << "] = " << vec(i) << ", index = " << indices(i) << endl;
    }
}