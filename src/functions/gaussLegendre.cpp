#include "../include/gaussLegendre.h"

using namespace std;
using namespace Eigen;

/**
 * @brief Calculate the Gauss-Legendre points and weights
 *
 * @param a lower bound of the integral
 * @param b upper bound of the integral
 * @param n number of points
 *
 * @return pair<VectorXd, VectorXd> nodes and weights
 */
pair<VectorXd, VectorXd> gaussLegendrePoints(double a, double b, size_t n) {
    VectorXd nodes(n), weights(n);
    if (n > 10) {
        gsl_integration_glfixed_table* table = gsl_integration_glfixed_table_alloc(n);
        double xi;
        double wi;
        for (int i = 0; i < n; ++i) {
            gsl_integration_glfixed_point(a, b, i, &xi, &wi, table);
            nodes(i) = xi;
            weights(i) = wi;
        }
        return { nodes, weights };
    }
    switch (n) {
    case 1:
        nodes << 0;
        weights << 2.0;
        break;
    case 2:
        nodes << 0.5773502691896257645091488, -0.5773502691896257645091488;
        weights << 1.0, 1.0;
        break;
    case 3:
        nodes << 0, 0.7745966692414833770358531, -0.7745966692414833770358531;
        weights << 0.8888888888888888888888889, 0.5555555555555555555555556, 0.5555555555555555555555556;
        break;
    case 4:
        nodes << 0.3399810435848562648026658, 0.8611363115940525752239465, -0.3399810435848562648026658, -0.8611363115940525752239465;
        weights << 0.6521451548625461426269361, 0.3478548451374538573730639, 0.6521451548625461426269361, 0.3478548451374538573730639;
        break;
    case 5:
        nodes << 0, 0.5384693101056830910363144, -0.5384693101056830910363144, 0.9061798459386639927976269, -0.9061798459386639927976269;
        weights << 0.5688888888888888888888889, 0.4786286704993664680412915, 0.4786286704993664680412915, 0.2369268850561890875142640, 0.2369268850561890875142640;
        break;
    case 6:
        nodes << 0.2386191860831969086305017, -0.2386191860831969086305017,
            0.6612093864662645136613996, -0.6612093864662645136613996,
            0.9324695142031520278123016, -0.9324695142031520278123016;
        weights << 0.4679139345726910473898703, 0.4679139345726910473898703,
            0.3607615730481386075698335, 0.3607615730481386075698335,
            0.1713244923791703450402961, 0.1713244923791703450402961;
        break;
    case 7:
        nodes << 0, -0.4058451513773971669066064, 0.4058451513773971669066064,
            0.7415311855993944398638648, -0.7415311855993944398638648,
            0.9491079123427585245261897, -0.9491079123427585245261897;
        weights << 0.4179591836734693877551020, 0.3818300505051189449503698,
            0.2797053914892766679014678, 0.1294849661688696932706114,
            0.3818300505051189449503698, 0.2797053914892766679014678,
            0.1294849661688696932706114;
        break;
    case 8:
        nodes << 0.1834346424956498049394761, -0.1834346424956498049394761,
            0.5255324099163289858177390, -0.5255324099163289858177390,
            0.7966664774136267395915539, -0.7966664774136267395915539,
            0.9602898564975362316835609, -0.9602898564975362316835609;
        weights << 0.3626837833783619829651504, 0.3626837833783619829651504,
            0.3137066458778872873379622, 0.3137066458778872873379622,
            0.2223810344533744705443560, 0.2223810344533744705443560,
            0.1012285362903762591525314, 0.1012285362903762591525314;
        break;
    case 9:
        nodes << 0, -0.3242534234038089290385380, 0.3242534234038089290385380,
            0.6133714327005903973087020, -0.6133714327005903973087020,
            0.8360311073266357942994298, -0.8360311073266357942994298,
            0.9681602395076260898355762, -0.9681602395076260898355762;
        weights << 0.3302393550012597631645251, 0.3302393550012597631645251,
            0.3123470770400028400686304, 0.3123470770400028400686304,
            0.2606106964029354623187429, 0.2606106964029354623187429,
            0.1806481606948574040584720, 0.1806481606948574040584720,
            0.0812743883615744119718922, 0.0812743883615744119718922;
        break;
    case 10:
        nodes << 0.1488743389816312108848260, -0.1488743389816312108848260,
            0.4333953941292471907992659, -0.4333953941292471907992659,
            0.6794095682990244062343274, -0.6794095682990244062343274,
            0.8650633666889845107320967, -0.8650633666889845107320967,
            0.9739065285171717200779640, -0.9739065285171717200779640;
        weights << 0.2955242247147528701738930, 0.2955242247147528701738930,
            0.2692667193099963550912269, 0.2692667193099963550912269,
            0.2190863625159820439955349, 0.2190863625159820439955349,
            0.1494513491505805931457763, 0.1494513491505805931457763,
            0.0666713443086881375935688, 0.0666713443086881375935688;
        break;
    default:
        break;
    }
    // normalize nodes and weights
    nodes.array() *= (b - a);
    nodes.array() += (a + b);
    nodes.array() /= 2.0;
    weights.array() *= (b - a) / 2.0;
    return { nodes, weights };
}

double gaussLegendre(std::function<double(double)> f, double a, double b, size_t n) {
    auto [nodes, weights] = gaussLegendrePoints(a, b, n); // Generalize
    // Calculate the integral using vectorized operations
    VectorXd fB = nodes.unaryExpr(f); // 
    double y = (fB.array() * weights.array()).sum();

    return y; // Return the integral result and the nodes
}

double gaussLegendre2D(std::function<double(double, double)> f, pair<double, double> x, pair<double, double> y, size_t n) {
    auto [xNodes, xWeights] = gaussLegendrePoints(x.first, x.second, n);
    auto [yNodes, yWeights] = gaussLegendrePoints(y.first, y.second, n);
    MatrixXd Values(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Values(i, j) = f(xNodes[i], yNodes[j]);
        }
    }

    // Calculate weights matrix
    MatrixXd Weights = xWeights * yWeights.transpose();
    Values.array() *= Weights.array(); // Element-wise multiplication    
    double result = Values.sum();

    return result; // Return the integral result
}

