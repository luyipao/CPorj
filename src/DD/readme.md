## 模型

$$ \begin{aligned}&n_{t}-(\mu En)_{x}=\tau\theta n_{xx},\\&\varphi_{xx}=\frac{e}{\varepsilon}(n-n_{d}),\end{aligned} $$

写成

$$ \sum_{1}^{N}\int_{I_{j}}n_{i}vdx+\sum_{1}^{N}(\int_{I_{j}}(\mu En)v_{x}dx+\int_{I_{j}}\tau\theta v_{x}dx)+\sum_{1}^{N}\left((\mu\widehat{En}+\tau\theta\widehat{n}_{x})[v]+\{v_{x}\}[\tau\theta n]\right)_{j+\frac{1}{2}},\quad n,v\in V_{h}^{k} $$

其中扩散通量
$$
\widehat{n_x} = \frac{\beta_0}{h} [n] + \set{n_x} + \beta_1 h [n_{xx}]
$$
通量
$$
\widehat{En} = \min\set{E,0}n^+ + \max \set{E, 0} n^-
$$

## 数值模拟

误差分析对于$\beta_0$敏感
