#include "../include/SolFunc2D.h"

#include <iomanip>

double normlizedP(int n, double x, double a, double b) {
    double map_x = (2 * x - a - b) / (b - a);
    double res = sqrt((2 * n + 1) / (b - a)) * gsl_sf_legendre_Pl(n, map_x);
    return res;
}

SolFunc2D::SolFunc2D(size_t maxDegree, double xa, double xb, double ya, double yb, size_t Nx, size_t Ny) {
    setNumBasis(maxDegree);
    this->grid.emplace_back(VectorXd::LinSpaced(Nx + 1, xa, xb));
    this->grid.emplace_back(VectorXd::LinSpaced(Ny + 1, ya, yb));
    this->cellNum = { Nx, Ny };
    this->coeff = vector<vector<VectorXd>>(Nx, vector<VectorXd>(Ny, VectorXd(this->basisNum)));
}

/**
 * @brief set information of the basis functions
 *
 * @param maxDegree basis function max degree
 * @param dim number of variables in the polynomial
 * @return int number of basis functions
 */
void SolFunc2D::setNumBasis(const size_t maxDegree) {
    this->maxDegree = maxDegree;
    this->basisNum = getBasisNum(2, maxDegree);
}
/**
 * @brief set the grid,
 *
 * @param x x direction grid
 * @param y y direction grid
 */
void SolFunc2D::setGrid(const vector<double>& x, const vector<double>& y) {
    this->cellNum = { x.size() - 1, y.size() - 1 };
    this->grid.emplace_back(VectorXd::Map(x.data(), x.size()));
    this->grid.emplace_back(VectorXd::Map(y.data(), y.size()));
}


void SolFunc2D::setCoeff(function<double(double, double)> f) {
    size_t n = 10;
    VectorXd xNodes(n), yNodes(n), xWeights(n), yWeights(n);
    MatrixXd fValues(n, n), Weight(n, n), basisValues(n, n);
    vector<int> varDegree(2);
    double res;
    double xa, xb, ya, yb;
    for (int i = 0; i < this->cellNum[0]; ++i) {
        for (int j = 0; j < this->cellNum[1]; ++j) {
            xa = this->grid[0](i);
            xb = this->grid[0](i + 1);
            ya = this->grid[1](j);
            yb = this->grid[1](j + 1);
            auto [xNodes, xWeights] = gaussLegendrePoints(xa, xb, n);
            auto [yNodes, yWeights] = gaussLegendrePoints(ya, yb, n);
            Weight = xWeights * yWeights.transpose();

            for (int k = 0; k < this->basisNum; ++k) {
                varDegree = getVarDegree(2, k);
                for (int ii = 0; ii < n; ++ii) {
                    for (int jj = 0; jj < n; ++jj) {
                        fValues(ii, jj) = f(xNodes(ii), yNodes(jj));
                        basisValues(ii, jj) = normlizedP(varDegree[0], xNodes(ii), xa, xb) *
                            normlizedP(varDegree[1], yNodes(jj), ya, yb);
                    }
                }
                res = (Weight.array() * fValues.array() * basisValues.array()).sum();
                this->coeff[i][j](k) = res;
                showProgressBar(this->cellNum[0] * this->cellNum[1] - 1, i * this->cellNum[1] + j);
            }
        }
    }
}
double SolFunc2D::operator()(vector<double> x) {

    for (int i = 0; i < 2; ++i) {
        if (x[i] < this->grid[i](0) && x[i] >= this->grid[i](this->cellNum[i])) {
            cerr << "x isn't in the grid" << endl;
        }
    }

    std::vector<int> position(2);
    for (int i = 0; i < 2; ++i) {
        auto l = std::lower_bound(this->grid[i].begin(), this->grid[i].end(), x[i]);

        position[i] = l - this->grid[i].begin();

        // 如果 position[i] 不为 0，减去 1
        if (position[i] > 0) {
            position[i] -= 1;
        }
    }

    double res = 0;
    for (int i = 0; i < this->basisNum; ++i) {
        vector<int> varDegree = getVarDegree(2, i);

        res += this->coeff[position[0]][position[1]][i] *
            normlizedP(varDegree[0], x[0], this->grid[0][position[0]], this->grid[0][position[0] + 1]) *
            normlizedP(varDegree[1], x[1], this->grid[1][position[1]], this->grid[1][position[1] + 1]);
    }
    return res;
}