#ifndef GAUSS_LEGENDRE_H
#define GAUSS_LEGENDRE_H
#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <functional>

using namespace std;
using namespace Eigen;

pair<VectorXd, VectorXd> gaussLegendrePoints(double a, double b, int n = 5);
double gaussLegendre(std::function<double(double)> f, double a, double b, int n = 5);
double gaussLegendre2D(std::function<double(double, double)> f, pair<double, double> x, pair<double, double> y, int n = 5);
#endif