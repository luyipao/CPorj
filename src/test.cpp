#include <iostream>
#include <Eigen/Dense>
#include "DDGIC.h"

using namespace Eigen;
using namespace std;

int main() {
    // Column-major 矩阵（默认）  
    DDGIC ddgic;
    ddgic.h = 1;
    ddgic.L();
    return 0;
}