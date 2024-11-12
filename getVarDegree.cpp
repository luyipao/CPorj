#include "getVarDegree.h"
#include "highDimFixedDegreePolyNum.h"
#include <vector>
#include <iostream>

using namespace std;

/**
 * @brief Get n-th corresponded Polynomial's degree as to every variable. n = 1,2,...
 *
 * @param dim polynomial's variable's num
 * @param n the n-th element of the polynomial, where n = 1, 2,...
 * @return vector<int> the degree of every variable, 0, 1, 2, ...
 */
vector<int> getVarDegree(int dim, int n) {
    int degree = 0;
    while (n > GhighDimFixedDegreePolyNum[dim][degree]) {
        n -= GhighDimFixedDegreePolyNum[dim][degree];
        ++degree;
    }
    vector<int> result(dim, 0);
    int i = 0;
    while (degree > 0) {
        while (n > GhighDimFixedDegreePolyNum[dim][degree - 1]) {
            n -= GhighDimFixedDegreePolyNum[dim][degree - 1];
            --dim;
            ++i;
        }
        while (degree > 0 && n <= GhighDimFixedDegreePolyNum[dim][degree - 1]) {
            degree--;
            ++result[i];
        }
    }
    return result;
}