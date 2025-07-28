// discretize.h

#ifndef DISCRETIZE_H
#define DISCRETIZE_H

#include <Eigen/Dense>
#include <stdexcept>
#include <limits>

/**
 * @brief 高效离散化（等距 edges 版本，C++ 0-based 索引）
 * @tparam DerivedX    任意 Eigen 1-D dense 实例 (VectorXd / ArrayXd / Block…)
 * @tparam DerivedEdge 同上，要求 monotone 且等距
 * @param  X           样本向量
 * @param  edges       区间边界（递增，等距），长度 N+1
 * @return VectorXi    与 X 同长度。
 *                     - 索引 0 到 N-1 
 *                     - 最后一个点 x == e_N 映射到 N-1
 *                     - 越界或 NaN 返回 -1 (或一个您选择的无效值)
 */
template <class DerivedX, class DerivedEdge>
inline Eigen::VectorXi discretize_uniform(const Eigen::DenseBase<DerivedX>& X,
    const Eigen::DenseBase<DerivedEdge>& edges) {
    using namespace Eigen;

    if (edges.size() < 2)
        throw std::invalid_argument("edges length must be ≥ 2");

    const double e0 = edges[0];
    const double step = edges[1] - edges[0];
    if (step <= 0.0)
        throw std::invalid_argument("edges must be strictly increasing");

    const int    nBins = static_cast<int>(edges.size()) - 1; // nBins 是区间的数量, 也是最大索引+1
    const double lastEdge = edges[edges.size() - 1];
    const double invStep = 1.0 / step;

    /* -------- 矢量化计算 -------- */
    // 核心计算：去掉 +1，现在直接映射到 0-based 索引
    // (x - e0) / step 的结果在 [0, 1) 内，floor 后是 0
    // (x - e0) / step 的结果在 [1, 2) 内，floor 后是 1
    ArrayXi idx =
        (((X.derived().array() - e0) * invStep).floor().template cast<int>());

    /* 右端闭：[e_{N-1}, e_N] —— x == lastEdge 映射到 nBins-1 */
    // 最后一个区间的索引是 nBins-1
    idx = (X.derived().array() == lastEdge).select(nBins - 1, idx);

    /* 越界 / NaN → -1 (一个典型的 C++ 无效索引值) */
    // 使用 -1 而不是 0，因为 0 现在是一个有效索引
    const int invalid_idx = -1; 
    auto bad = (X.derived().array() <  e0) ||
               (X.derived().array() >  lastEdge) ||
               (X.derived().array().isNaN());
    idx = bad.select(invalid_idx, idx);

    /* clip，防止数值误差导致 nBins (越界) */
    // 有效索引范围是 [0, nBins-1]
    idx = idx.min(nBins - 1);

    return idx.matrix();
}

#endif // DISCRETIZE_H