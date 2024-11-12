import numpy as np
import time

# 定义两个向量
N = int(1e7)  # 将 N 转换为整数

a = np.random.rand(N)
b = np.random.rand(N)  

# 使用 NumPy 的 dot 函数计算点积
start_time = time.time()
dot_product = np.dot(a, b)
end_time = time.time()
numpy_time = (end_time - start_time)  # 转换为毫秒

# 使用循环计算点积
result = 0
start_time = time.time()
for i in range(N):
    result += a[i] * b[i]
end_time = time.time()
loop_time = (end_time - start_time) # 转换为毫秒
del a
del b
# 输出结果
print("NumPy 点积:", dot_product)
print("循环计算点积:", result)
print("NumPy 计算时间: {:.6f} s".format(numpy_time))  # 输出毫秒
print("循环计算时间: {:.6f} s".format(loop_time))  # 输出毫秒
