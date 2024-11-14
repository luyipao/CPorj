#include "include/LegendrePoly.h"
#include "include/getVarDegree.h"
#include <cmath>

int main() {
    vector<vector<double>> mesh = { {-1, 0,1},{-1, 0,1} };
    LegendrePoly P(2, 3, mesh);
    P.getCoeff();
    
    cout << P({ 0,0 }) << ' ' << P({ 1,1 }) << ' ' << P({ -1,-1 }) << endl;
    cout << P({0.5,0.5}) << ' ' << P({ 0.5,-0.5 }) << ' ' << P({ -0.5,0.5 }) << ' ' << P({ -0.5,-0.5 }) << endl;
    

    return 0;
}