#ifndef GETBASISNUM_H
#define GETBASISNUM_H
/**
 * @brief get the number of basis functions
 *
 * @param maxDegree basis function max degree
 * @param dim number of variables in the polynomial
 * @return int number of basis functions
 */
long long int getBasisNum(const int dim, const int maxDegree);

#endif