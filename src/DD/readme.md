## 模型

$$
\begin{aligned}

&n_{t}-(\mu En)_{x}=\tau\theta n_{xx},\\
&\varphi_{xx}=\frac{e}{\varepsilon}(n-n_{d}), \\
& E = -\phi_x
\end{aligned}
$$

其中未知量是$n,\phi,E$，其余为常量（$\mathcal{k}$表示玻尔兹曼常量）。假设$n$光滑，且具有周期边界条件。
$$
\mathcal{k} = 0.138e-4,\varepsilon = 11.7\times 8.85418,e = 0.1602,m = 0.26\times0.9109e-31,\mu = 0.75,T_0 = 300, \\
\tau = \frac{m\mu}{e}, \theta = \frac{\mathcal{k}}{m}T_0
$$
初值条件，使用$\mathrm{S}_4(x)=70x^9-315x^8+540x^7-420x^6+126x^5$过渡
$$
n_d = 
\begin{cases}
5e5,x \in [0,0.1]\\
5e5\left(1-S_4(\frac{x-0.1}{0.05})\right) +2e3 S_4(\frac{x-0.1}{0.05}),x\in(0.1,0.15) \\
2e3, x \in [0.15, 0.45] \\
2e3\left(1-S_4(\frac{x-0.45}{0.05})+5e5S_4(\frac{x-0.45}{0.05}) \right),x\in(0.45,0.5) \\
5e5,x \in[0.5,0.6] 

\end{cases}
$$

E采用连续性方法
$$
E^h=\int_0^x-\frac e\varepsilon(n^h-n_d)dx+E_0-v_\mathrm{bias},E_0\:=\:E^h(0)\:=\:\int_0^1(\int_0^x\frac e\varepsilon(n^h-n_d)ds)dx.
$$
数值格式
$$
\int_{I_j} (n^h)_t v - (\mu \widehat{E^hn^h} + \tau\theta \widehat{n^h_x})v\Big|_{I_j} + \int( \mu E^h n^h + \tau\theta (n^h)_x )v_x +\frac{1}{2}[b(n^h)](v_x)_{j+\frac{1}{2}}^- + \frac{1}{2}[b(u)](v_x)_{j-\frac{1}{2}}^+ = 0
$$


其中扩散通量
$$
\widehat{n_x} = \frac{\beta_0}{h} [n] + \set{n_x} + \beta_1 h [n_{xx}]
$$
通量
$$
\widehat{En} = \min\set{E,0}n^- + \max \set{E, 0} n^+
$$

$$
b[n] = \tau \theta n
$$

## 基函数

定义区间$I_{j}=\left[ x_{j-\frac{1}{2}},x_{j+\frac{1}{2}} \right]$上的基函数
$$
\phi_{m} = \left( \frac{\left( x-x_{j-\frac{1}{2}} \right)}{\Delta x} \right)^{m}, m = 0,1,2,\dots
基函数的导数
$$

$$
\phi_{0}' = 0,\quad \phi'_{m} = \frac{m}{\Delta x} \phi_{m-1},m\geq 1
$$
则我们的目标函数可以定义为
$$
u = \sum_{m=0}^k c_{m}^j \phi ^{m}, \quad \phi = \frac{ x-x_{j-\frac{1}{2}} }{\Delta x}
$$
## 数值模拟

目前没有
