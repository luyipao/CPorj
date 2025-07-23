#ifndef PHYCONST_H
#define PHYCONST_H
class PhyConst {
public:
    static constexpr double mu = 0.75;
    static constexpr double T_0 = 300;
    static constexpr double m = 0.26 * 0.9109e-31; 
    static constexpr double e = 0.1602;           
    static constexpr double k = 0.138e-4;        
    static constexpr double epsilon = 11.7 * 8.85418; 

    static constexpr double tau = m * mu / e;
    static constexpr double theta = k * T_0 / m;
};

#endif