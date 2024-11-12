#include <iostream>
#include <Eigen/Dense>
#include <chrono>
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

int s() {
    int N = 1e7;

    VectorXd a(N);
    VectorXd b(N);
    a.setRandom();
    b.setRandom();
    double result;

    auto start = chrono::high_resolution_clock::now();
    result = a.dot(b);
    auto end = chrono::high_resolution_clock::now();


    chrono::duration<double> elapsed = end - start;
    cout << elapsed.count() << " seconds" << endl;
    cout << result << endl;
    result = 0;
    start = chrono::high_resolution_clock::now();
    for (int i = 0; i < N; i++) {
        result += a(i) * b(i);
    }
    end = chrono::high_resolution_clock::now();
    elapsed = end - start;
    cout << elapsed.count() << " seconds" << endl;
    cout << result << endl;
    return 0;
}