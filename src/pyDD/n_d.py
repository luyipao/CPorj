import numpy as np
import matplotlib.pyplot as plt

def smooth_step(x):
    """
    平滑步进函数
    支持标量和数组输入
    """
    if np.isscalar(x):
        return 70 * x**9 - 315 * x**8 + 540 * x**7 - 420 * x**6 + 126 * x**5
    else:
        # 数组版本
        x = np.asarray(x)
        return 70 * x**9 - 315 * x**8 + 540 * x**7 - 420 * x**6 + 126 * x**5

def n_d(X):
    """
    计算n_d函数值
    
    参数:
    X: 输入点数组或标量
    
    返回:
    对应的n_d值
    """
    # 处理标量输入
    if np.isscalar(X):
        X = np.array([X])
        return_scalar = True
    else:
        X = np.asarray(X)
        return_scalar = False
    
    Y = np.zeros(X.shape)
    
    for i, x in enumerate(X):
        if x < 0 or x > 0.6:
            raise ValueError(f"X value out of bounds: {x}")
        elif x >= 0 and x <= 0.1:
            Y[i] = 5e5
        elif x > 0.1 and x < 0.15:
            x_norm = (x - 0.1) / 0.05
            Y[i] = 5e5 * (1 - smooth_step(x_norm)) + 2e3 * smooth_step(x_norm)
        elif x >= 0.15 and x <= 0.45:
            Y[i] = 2e3
        elif x > 0.45 and x < 0.5:
            x_norm = (x - 0.45) / 0.05
            Y[i] = 2e3 * (1 - smooth_step(x_norm)) + 5e5 * smooth_step(x_norm)
        elif x >= 0.5 and x <= 0.6:
            Y[i] = 5e5
    
    return Y[0] if return_scalar else Y

def main():
    """
    主函数 - 生成数据并保存到文件
    """
    X = np.linspace(0, 0.6, 10000)
    Y = n_d(X)
    
    # 保存数据到文件
    with open('data.txt', 'w') as file:
        for i in range(len(X)):
            file.write(f"{X[i]} {Y[i]}\n")
    print("数据已保存到 data.txt")
    
    # 可选：绘制图形
    plt.figure(figsize=(10, 6))
    plt.plot(X, Y, 'b-', linewidth=2)
    plt.xlabel('x')
    plt.ylabel('n_d(x)')
    plt.title('n_d Function')
    plt.grid(True)
    plt.yscale('log')  # 由于数值范围很大，使用对数坐标
    plt.show()

if __name__ == "__main__":
    main()