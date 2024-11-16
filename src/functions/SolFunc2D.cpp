#include "../include/SolFunc2D.h"
#include "../include/getBasisNum.h"
#include <iomanip>


/**
 * @brief set information of the basis functions
 *
 * @param maxDegree basis function max degree
 * @param dim number of variables in the polynomial
 * @return int number of basis functions
 */
void SolFunc2D::setNumBasis(const int maxDegree) {
    this->maxDegree = maxDegree;
    this->basisNum = getBasisNum(2, maxDegree);
}
/**
 * @brief set the mesh,
 *
 * @param x x direction mesh
 * @param y y direction mesh
 */
void SolFunc2D::setGrid(const vector<double>& x, const vector<double>& y) {
    this->cellNum = { x.size() - 1, y.size() - 1 };
    this->grid.emplace_back(VectorXd::Map(x.data(), x.size()));
    this->grid.emplace_back(VectorXd::Map(y.data(), y.size()));
}


