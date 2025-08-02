import numpy as np
import datetime
import math

import n_d
import DDGIC4DD 

def main():
    start = datetime.datetime.now()
    start_time_str = start.strftime('%c') # 格式化开始时间

    k = 2
    Xa = 0.00
    Xb = 0.60

    # double beta0 = (k + 1) * (k + 1);
    # double beta1 = 1 / (k + 1.0) / (2 * k);
    beta0 = 32.0
    beta1 = 1.0 / 12.0
    T = 1.0
    
    if k == 0:
        CFL = 0.1
    elif k == 1:
        CFL = 0.1
    else:
        CFL = 0.1

    # 构建初始函数系数
    Ns = [80, 160, 320]  # 网格数
    h_values = []
    # C++ 代码中 L2_errors 和 LInf_errors 的最后一个元素不会被赋值，
    # np.zeros 的行为与此一致。
    L2_errors = np.zeros(len(Ns))
    LInf_errors = np.zeros(len(Ns))

    for N in Ns:
        h = (Xb - Xa) / N
        mesh = np.linspace(Xa, Xb, N + 1)  # Mesh points
        n_C = np.zeros((k + 1, N))        # Coefficient matrix

        for j in range(N):
            # 构建质量矩阵 M
            M = np.zeros((k + 1, k + 1))
            for i in range(k + 1):
                for l in range(k + 1):
                    M[i, l] = h / (i + l + 1)

            # 使用足够多的积分点以精确投影初始条件
            nodes, weights = gaussLegendrePoints(mesh[j], mesh[j + 1], n=k + 5)

            b = np.zeros(k + 1)
            for i in range(k + 1):
                # 将初始条件函数 n_d 投影到基函数上
                b[i] = np.sum(weights * n_d(nodes) * (((nodes - mesh[j]) / h) ** i))

            # 求解线性系统 M * c = b
            n_C[:, j] = np.linalg.solve(M, b)
        
        ddgic = DDGIC4DD(N, Xa, Xb, k, beta0, beta1, T, CFL)
        # ddgic.dt = 0.0000225 # 可以取消注释以覆盖 dt
        ddgic.setn(n_C)
        ddgic.setn_d(n_C) # 原始代码也将 n_C 设置给了 n_d_C

        # 执行时间步进
        ddgic.RKDDG()

        # 在一个精细网格上评估最终解，用于误差分析
        X_fine = np.linspace(Xa, Xb, 10001)
        Y_fine = ddgic.getn(X_fine)
        h_values.append(Y_fine)

    # --- 误差分析 ---
    filename = "ErrorAnalysis.txt"
    try:
        # 使用 'a' 模式追加内容到文件
        with open(filename, 'a', encoding='utf-8') as result_file:
            end = datetime.datetime.now()
            
            # --- 2. 获取并写入当前时间戳 ---
            result_file.write("======================================================================\n")
            # C++ 的 ctime() 会自带换行符，这里我们也添加 '\n'
            result_file.write(f"开始于: {start_time_str}\n结束于：{end.strftime('%c')}\n")
            result_file.write("======================================================================\n")

            result_file.write(f"k = {k}, beta0 = {beta0}, beta1 = {beta1:.6f}, T = {T}, CFL = {CFL}\n\n")

            for i in range(len(Ns) - 1):
                # 解都在同一个精细网格上，所以可以直接相减
                # 这计算的是 N_i 和 N_{i+1} 两种网格密度下的解的差
                error_vec = h_values[i] - h_values[i + 1]
                h_i = (Xb - Xa) / Ns[i]
                
                # L2 范数计算，与 C++ 代码完全对应
                # error_vec.size() 对应 len(error_vec)
                L2_errors[i] = np.linalg.norm(error_vec) * math.sqrt((Xb - Xa) / len(error_vec))
                
                # 无穷范数计算
                LInf_errors[i] = np.linalg.norm(error_vec, ord=np.inf)

                result_file.write(f"N: {Ns[i]} (h={h_i:.6f}), L2 Error: {L2_errors[i]:.6e}, LInf Error: {LInf_errors[i]:.6e}\n")

            # --- 计算收敛阶 ---
            result_file.write("\n--- 收敛阶 ---\n")
            # 注意循环边界，避免索引越界
            for i in range(len(Ns) - 2):
                # 收敛阶 = log(E_h / E_{h/2}) / log(2)
                l2_order = math.log(L2_errors[i] / L2_errors[i+1]) / math.log(2.0)
                l_inf_order = math.log(LInf_errors[i] / LInf_errors[i+1]) / math.log(2.0)

                result_file.write(f"Order from N={Ns[i]} to N={Ns[i+1]}: L2 Order = {l2_order:.6f}, L_inf Order = {l_inf_order:.6f}\n")
            
            result_file.write("\n\n")
        
        print(f"计算完成，结果已追加到 {filename}")

    except IOError as e:
        print(f"错误: 无法打开结果文件 {filename}: {e}")

if __name__ == "__main__":
    main()