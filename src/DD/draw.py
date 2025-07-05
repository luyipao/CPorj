import numpy as np
import matplotlib.pyplot as plt

# 加载数据
data = np.loadtxt('data.txt')
x = data[:, 0]
y = data[:, 1]

# 创建图形和轴
fig, ax = plt.subplots(figsize=(10, 6), dpi=120)  # 提高DPI获得更精细图像

# 绘制函数曲线（更精细的参数）
line, = ax.plot(x, y, 
                linewidth=1.5,  # 增加线宽
                color='royalblue',
                antialiased=True,  # 启用抗锯齿
                zorder=1, 
                label='n_d(x)')

# 设置标题和标签
ax.set_title('Function: n_d(x)', fontsize=14)
ax.set_xlabel('x', fontsize=12)
ax.set_ylabel('n_d(x)', fontsize=12)

# 设置更专业的网格
ax.grid(True, linestyle='--', linewidth=0.7, alpha=0.7, which='both')

# 添加图例
ax.legend(loc='best', framealpha=0.8)

# 添加交互式点选择功能
def on_click(event):
    # 检查是否在坐标轴内点击
    if event.inaxes == ax:
        # 找到最近的数据点
        x_click = event.xdata
        idx = np.argmin(np.abs(x - x_click))
        x_point = x[idx]
        y_point = y[idx]
        
        # 清除旧标记
        if hasattr(on_click, 'marker'):
            on_click.marker.remove()
            on_click.annotation.remove()
        
        # 绘制新标记
        on_click.marker = ax.scatter(x_point, y_point, 
                                     color='red', 
                                     s=60, 
                                     edgecolor='black',
                                     zorder=3,
                                     label='Selected Point')
        
        # 添加坐标文本
        on_click.annotation = ax.annotate(f'({x_point:.3f}, {y_point:.3f})',
                                          (x_point, y_point),
                                          xytext=(10, -15),
                                          textcoords='offset points',
                                          fontsize=10,
                                          bbox=dict(boxstyle="round,pad=0.3", 
                                                    fc="white", 
                                                    ec="gray", 
                                                    alpha=0.8),
                                          arrowprops=dict(arrowstyle="->", 
                                                          connectionstyle="arc3"))
        
        # 刷新图形
        fig.canvas.draw_idle()
        print(f"Selected point: x={x_point:.6f}, y={y_point:.6f}")

# 连接点击事件
fig.canvas.mpl_connect('button_press_event', on_click)
plt.tight_layout()  # 优化布局

# 保存高分辨率图像
plt.savefig('function_plot.png', dpi=300, bbox_inches='tight')

# 显示图形
plt.show()