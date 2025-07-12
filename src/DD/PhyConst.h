#ifndef PHYCONST_H
#define PHYCONST_H

class PhyConst {
public:
    static constexpr double mu = 0.75;
    static constexpr double T_0 = 300;
    static constexpr double m = 0.26 * 0.9109e-31;   // 修正电子质量
    static constexpr double e = 0.1602;           // 修正基本电荷
    static constexpr double k = 0.138 * 1e-4;        // 玻尔兹曼常数
    static constexpr double epsilon = 11.7 * 8.85418; // 真空介电常数

    // 计算相关常量
    static constexpr double tau = m * mu / e;
    static constexpr double theta = k * T_0 / m;
};

#endif // PHYCONST_H