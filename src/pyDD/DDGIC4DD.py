import numpy as np
from scipy.sparse import csc_matrix, diags
from scipy.sparse.linalg import spsolve
from scipy.linalg import solve
import matplotlib.pyplot as plt
from scipy import integrate
from typing import Tuple, Dict, List
import warnings

class PhyConst:
    """物理常数类"""
    mu = 0.75
    T_0 = 300
    m = 0.26 * 0.9109e-31
    e = 0.1602
    k = 0.138e-4
    epsilon = 11.7 * 8.85418
    
    tau = m * mu / e
    theta = k * T_0 / m

def gaussLegendrePoints(a: float, b: float, n: int = 5) -> Tuple[np.ndarray, np.ndarray]:
    """
    计算区间[a,b]上的高斯-勒让德积分点和权重
    """
    nodes, weights = np.polynomial.legendre.leggauss(n)
    # 转换到区间[a,b]
    nodes = 0.5 * (b - a) * nodes + 0.5 * (a + b)
    weights = 0.5 * (b - a) * weights
    return nodes, weights

def discretize_uniform(X: np.ndarray, mesh: np.ndarray) -> np.ndarray:
    """
    将点X离散化到均匀网格上，返回每个点所在的单元格索引
    """
    h = mesh[1] - mesh[0]
    Xa = mesh[0]
    indices = np.floor((X - Xa) / h).astype(int)
    # 处理边界情况
    indices = np.clip(indices, 0, len(mesh) - 2)
    # 处理超出边界的点
    out_of_bounds = (X < mesh[0]) | (X > mesh[-1])
    indices[out_of_bounds] = -1
    return indices

class DDGIC4DD:
    def __init__(self, N: int, Xa: float, Xb: float, k: int, 
                 beta0: float, beta1: float, T: float, CFL: float):
        """
        初始化DDGIC4DD类
        
        参数:
        N: 网格单元数
        Xa: 左边界
        Xb: 右边界
        k: 多项式最大阶数
        beta0, beta1: 参数
        T: 总时间
        CFL: CFL条件
        """
        self.N = N
        self.Xa = Xa
        self.Xb = Xb
        self.k = k
        self.beta0 = beta0
        self.beta1 = beta1
        self.T = T
        self.CFL = CFL
        
        # 计算网格参数
        self.h = (Xb - Xa) / N
        self.dt = CFL * self.h * self.h  # CFL条件
        self.mesh = np.linspace(Xa, Xb, N + 1)  # 网格点
        
        # 预计算高斯积分点
        self.gaussPoints: Dict[int, Tuple[np.ndarray, np.ndarray]] = {}
        for j in range(N):
            self.gaussPoints[j] = gaussLegendrePoints(self.mesh[j], self.mesh[j + 1])
        
        # 构建质量矩阵A
        self._build_mass_matrix()
        
        # L()函数中用到的常量
        self.LL, self.LR, self.RL, self.RR = 0, 1, 2, 3
        
        # 用于存储导数计算的帮助向量
        self.uxHelp = np.arange(k + 1, dtype=float)
        self.uxxHelp = self.uxHelp * (self.uxHelp - 1)
        
        # 初始化系数矩阵
        self.C = np.zeros((k + 1, N))
        self.n_d_C = np.zeros((k + 1, N))
    
    def _build_mass_matrix(self):
        """构建质量矩阵"""
        # 构建块对角质量矩阵
        size = (self.k + 1) * self.N
        data = []
        row_indices = []
        col_indices = []
        
        for j in range(self.N):
            block_offset = j * (self.k + 1)
            for l in range(self.k + 1):
                for m in range(self.k + 1):
                    value = self.h / (m + l + 1.0)
                    global_row = block_offset + l
                    global_col = block_offset + m
                    data.append(value)
                    row_indices.append(global_row)
                    col_indices.append(global_col)
        
        self.A_sparse = csc_matrix((data, (row_indices, col_indices)), shape=(size, size))
    
    def setn(self, coeff: np.ndarray):
        """
        设置n的系数矩阵
        
        参数:
        coeff: 系数矩阵，每列代表一个单元的系数
        """
        if coeff.shape != (self.k + 1, self.N):
            raise ValueError("Coefficient matrix size mismatch.")
        self.C = coeff.copy()
    
    def setn_d(self, coeff: np.ndarray):
        """
        设置n_d的系数矩阵
        
        参数:
        coeff: 系数矩阵，每列代表一个单元的系数
        """
        if coeff.shape != (self.k + 1, self.N):
            raise ValueError("Coefficient matrix size mismatch.")
        self.n_d_C = coeff.copy()
    
    def getPriE_x(self, X: np.ndarray) -> np.ndarray:
        """
        计算E_x的原函数（从Xa到X的定积分）
        
        参数:
        X: 要求积分的点
        
        返回:
        在每个点X(i)处的定积分
        """
        if np.isscalar(X):
            X = np.array([X])
        
        Coeff = self.C - self.n_d_C  # 被积函数的多项式系数
        Y = np.zeros(X.shape)
        
        # 预计算系数矩阵的累积和
        CumSumCoeff = np.cumsum(Coeff, axis=1)
        
        # 获取所有点所在的单元格索引
        bins = discretize_uniform(X, self.mesh)
        
        for i, xi in enumerate(X):
            if xi == self.Xa:
                Y[i] = 0.0
                continue
                
            bin_idx = bins[i]
            if bin_idx < 0:
                Y[i] = np.nan
                continue
            
            u = (xi - self.mesh[bin_idx]) / self.h
            integral_at_xi = 0.0
            
            # 对每个多项式阶数m进行计算
            for m in range(self.k + 1):
                m_plus_1 = m + 1
                full_cells_integral_contrib = 0.0
                
                if bin_idx > 0:
                    full_cells_integral_contrib = CumSumCoeff[m, bin_idx - 1] * self.h / m_plus_1
                
                partial_cell_integral_contrib = Coeff[m, bin_idx] * self.h * (u ** m_plus_1) / m_plus_1
                integral_at_xi += full_cells_integral_contrib + partial_cell_integral_contrib
            
            Y[i] = integral_at_xi
        
        return Y if len(Y) > 1 else Y[0]
    
    def discretize(self, x: float) -> int:
        """将点x离散化到网格上"""
        if x < self.Xa or x > self.Xb:
            raise ValueError("x is out of bounds.")
        if x == self.Xb:
            return self.N - 1
        return int((x - self.Xa) / self.h)
    
    def getE(self, X: np.ndarray) -> np.ndarray:
        """
        获取电场E在给定点X处的值
        """
        nodes, weights = gaussLegendrePoints(0, 0.6)
        E0 = np.dot(self.getPriE_x(nodes), weights)
        E0 *= PhyConst.e / PhyConst.epsilon
        
        y = -PhyConst.e / PhyConst.epsilon * self.getPriE_x(X) + E0 - 1.5
        return y
    
    def getn(self, x) -> np.ndarray:
        """
        获取n在给定点处的值
        """
        if np.isscalar(x):
            # 标量版本
            idx = self.discretize(x)
            Cellx = (x - self.mesh[idx]) / self.h
            exp_terms = np.array([Cellx**i for i in range(self.k + 1)])
            return np.dot(self.C[:, idx], exp_terms)
        else:
            # 向量版本
            X = np.asarray(x)
            M = X.size
            idx = discretize_uniform(X, self.mesh)
            
            # 计算每个单元内的归一化坐标
            leftEdges = self.Xa + idx * self.h
            Cellx = (X - leftEdges) / self.h
            
            # 使用Horner方法进行多项式求值
            Y = np.zeros(M)
            for j in range(self.k, -1, -1):
                # 为每个点的单元索引收集第j个系数
                coeff_j = np.array([self.C[j, cell] if cell >= 0 else 0 for cell in idx])
                Y = Y * Cellx + coeff_j
            
            return Y
    
    def getn_x(self, x) -> np.ndarray:
        """
        获取n的一阶导数在给定点处的值
        """
        if np.isscalar(x):
            # 标量版本
            idx = self.discretize(x)
            Cellx = (x - self.mesh[idx]) / self.h
            y = 0.0
            tempx = 1.0
            for i in range(1, self.k + 1):
                y += self.C[i, idx] * i * tempx
                tempx *= Cellx
            return y / self.h
        else:
            # 向量版本
            X = np.asarray(x)
            result = np.zeros(X.shape)
            for i, xi in enumerate(X):
                result[i] = self.getn_x(xi)  # 递归调用标量版本
            return result
    
    def phi_x(self, X: np.ndarray, k: int) -> np.ndarray:
        """
        计算基函数的导数
        """
        X = np.asarray(X)
        Y = np.zeros(X.shape)
        
        for i, xi in enumerate(X):
            idx = self.discretize(xi)
            Cellx = (xi - self.mesh[idx]) / self.h
            Y[i] = k / self.h * (Cellx ** (k - 1)) if k > 0 else 0
        
        return Y
    
    def L(self) -> np.ndarray:
        """
        计算L操作符
        """
        # 初始化u, ux, uxx
        u = np.zeros((4, self.N))
        ux = np.zeros((4, self.N))
        uxx = np.zeros((4, self.N))
        
        # 计算u值
        if self.N > 1:
            u[self.LL, 1:] = np.sum(self.C[:, :-1], axis=0)  # LL
            u[self.RR, :-1] = self.C[0, 1:]  # RR
        u[self.LR, :] = self.C[0, :]
        u[self.RL, :] = np.sum(self.C, axis=0)
        
        # 计算ux值
        if self.k >= 1:
            C_dot_uxHelp = np.dot(self.C.T, self.uxHelp)
            if self.N > 1:
                ux[self.LL, 1:] = C_dot_uxHelp[:-1]  # LL
                ux[self.RR, :-1] = self.C[1, 1:]  # RR
            ux[self.LR, :] = self.C[1, :]
            ux[self.RL, :] = C_dot_uxHelp
            ux = ux / self.h
        
        # 计算uxx值
        if self.k >= 2:
            C_dot_uxxHelp = np.dot(self.C.T, self.uxxHelp)
            if self.N > 1:
                uxx[self.LL, 1:] = C_dot_uxxHelp[:-1]  # LL
                uxx[self.RR, :-1] = 2 * self.C[2, 1:]
            uxx[self.LR, :] = 2 * self.C[2, :]
            uxx[self.RL, :] = C_dot_uxxHelp
            uxx = uxx / (self.h * self.h)
        
        # 周期边界条件
        if self.k >= 0:
            u[self.LL, 0] = u[self.RL, self.N - 1]
            u[self.RR, self.N - 1] = u[self.LR, 0]
        if self.k >= 1:
            ux[self.LL, 0] = ux[self.RL, self.N - 1]
            ux[self.RR, self.N - 1] = ux[self.LR, 0]
        if self.k >= 2:
            uxx[self.LL, 0] = uxx[self.RL, self.N - 1]
            uxx[self.RR, self.N - 1] = uxx[self.LR, 0]
        
        # 预先计算所有需要的E值
        num_mesh_points = len(self.mesh)
        num_gauss_nodes = len(self.gaussPoints[0][0])
        total_points = num_mesh_points + self.N * num_gauss_nodes
        
        combined_points = np.zeros(total_points)
        combined_points[:num_mesh_points] = self.mesh
        
        for j in range(self.N):
            nodes, weights = self.gaussPoints[j]
            start_idx = num_mesh_points + j * num_gauss_nodes
            combined_points[start_idx:start_idx + num_gauss_nodes] = nodes
        
        E_all = self.getE(combined_points)
        E = E_all[:num_mesh_points]
        E_gauss = []
        for j in range(self.N):
            start_idx = num_mesh_points + j * num_gauss_nodes
            E_gauss.append(E_all[start_idx:start_idx + num_gauss_nodes])
        
        # 预计算phi_x_gauss
        phi_x_gauss = []
        for j in range(self.N):
            nodes, _ = self.gaussPoints[j]
            phi_x_cell = []
            for i in range(self.k + 1):
                phi_x_cell.append(self.phi_x(nodes, i))
            phi_x_gauss.append(phi_x_cell)
        
        # 计算通量项B
        B = np.zeros((self.k + 1, self.N))
        
        tempB1 = PhyConst.mu * (
            np.minimum(E[1:self.N+1], 0.0) * u[self.RL, :] +
            np.maximum(E[1:self.N+1], 0.0) * u[self.RR, :]
        )
        
        tempB2 = self.beta0 / self.h * (u[self.RR, :] - u[self.RL, :])
        tempB3 = 0.5 * (ux[self.RR, :] + ux[self.RL, :])
        tempB4 = self.beta1 * self.h * (uxx[self.RR, :] - uxx[self.RL, :])
        
        for l in range(self.k + 1):
            B[l, :] = tempB1 + PhyConst.tau * PhyConst.theta * (tempB2 + tempB3 + tempB4)
            
            if l == 0:
                tempB5 = PhyConst.mu * (
                    np.minimum(E[:self.N], 0.0) * u[self.LL, :] +
                    np.maximum(E[:self.N], 0.0) * u[self.LR, :]
                )
                
                tempB6 = self.beta0 * (u[self.LR, :] - u[self.LL, :]) / self.h
                tempB7 = 0.5 * (ux[self.LR, :] + ux[self.LL, :])
                tempB8 = self.beta1 * self.h * (uxx[self.LR, :] - uxx[self.LL, :])
                
                B[l, :] -= tempB5 + PhyConst.tau * PhyConst.theta * (tempB6 + tempB7 + tempB8)
        
        # 计算积分项D
        D = np.zeros((self.k + 1, self.N))
        for j in range(self.N):
            nodes, weights = self.gaussPoints[j]
            for l in range(self.k + 1):
                D[l, j] = PhyConst.mu * np.sum(weights * E_gauss[j] * self.getn(nodes) * phi_x_gauss[j][l])
                D[l, j] += PhyConst.tau * PhyConst.theta * np.sum(weights * self.getn_x(nodes) * phi_x_gauss[j][l])
        
        # 计算界面修正项E
        E_mat = np.zeros((self.k + 1, self.N))
        if self.k >= 1:
            for l in range(1, self.k + 1):
                for j in range(self.N):
                    term1 = (u[self.RR, j] - u[self.RL, j]) * l / self.h
                    term2 = (u[self.LR, j] - u[self.LL, j]) * l / self.h * (1 if l == 1 else 0)
                    E_mat[l, j] = PhyConst.tau * PhyConst.theta * 0.5 * (term1 + term2)
        
        # 构建右端向量
        b_mat = B - D - E_mat
        b = b_mat.flatten('F')  # 按列展开，对应Eigen的默认存储顺序
        
        # 求解线性系统
        res_vec = spsolve(self.A_sparse, b)
        
        # 重塑结果为矩阵
        res = res_vec.reshape((self.k + 1, self.N), order='F')
        
        return res
    
    def RK(self):
        """
        Runge-Kutta时间步进
        """
        C_pre = self.C.copy()
        
        # k1步
        k1 = self.L()
        k1 = C_pre + self.dt * k1
        self.setn(k1)
        
        # k2步
        k2 = self.L()
        k2 = 3.0/4.0 * C_pre + 1.0/4.0 * k1 + 1.0/4.0 * self.dt * k2
        self.setn(k2)
        
        # k3步
        k3 = self.L()
        k3 = 1.0/3.0 * C_pre + 2.0/3.0 * k2 + 2.0/3.0 * self.dt * k3
        self.setn(k3)
    
    def RKDDG(self):
        """
        运行完整的Runge-Kutta DDG求解
        """
        TT = 0
        Sols = []
        X = np.linspace(self.Xa, self.Xb, 10001)
        
        while TT < self.T:
            # Y = self.getn(X)
            # Sols.append(Y)
            self.RK()
            TT += self.dt
    
    def drawn(self):
        """
        绘制结果并保存数据
        """
        num_points = 1000
        x = np.linspace(0, 0.6, num_points)
        y = self.getn(x)
        
        # 保存数据到文件
        with open('data.txt', 'w') as file:
            for i in range(len(x)):
                file.write(f"{x[i]} {y[i]}\n")
        print("数据已保存到 data.txt")
        
        # 可选：绘制图形
        plt.figure(figsize=(10, 6))
        plt.plot(x, y, 'b-', linewidth=2)
        plt.xlabel('x')
        plt.ylabel('n(x)')
        plt.title('Solution n(x)')
        plt.grid(True)
        plt.show()

# 使用示例
if __name__ == "__main__":
    # 创建求解器实例
    solver = DDGIC4DD(N=100, Xa=0.0, Xb=0.6, k=2, 
                     beta0=1.0, beta1=0.1, T=1.0, CFL=0.1)
    
    # 设置初始条件（这里需要根据具体问题设置）
    initial_coeff = np.random.rand(solver.k + 1, solver.N) * 0.1
    solver.setn(initial_coeff)
    solver.setn_d(np.zeros((solver.k + 1, solver.N)))
    
    # 运行求解
    # solver.RKDDG()
    
    # 绘制结果
    solver.drawn()