#include "../include/getBasisNum.h"
long long int getBasisNum(const int dim, const int maxDegree) {
    long long int a = 1;
    long long int b = 1;
    for (int i = 1; i <= maxDegree; i++) {
        a *= (dim + maxDegree);
        b *= i;
    }
    return a / b;
}