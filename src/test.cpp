#include <iostream>
#include <Eigen/Dense>
#include "DDGIC.h"

using namespace Eigen;
using namespace std;

int main() {
    // Column-major 矩阵（默认）  
    int N = 10;
    VectorXd mesh = VectorXd::LinSpaced(N + 1, 0, 1); // Mesh points
    double h = (mesh(N) - mesh(0)) / N; // 网格大小
    double x = 1;
    double Xa = 0;
    cout << mesh.size() << endl;
    return 0;
}