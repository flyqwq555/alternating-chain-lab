# alternating-chain-lab

交替高斯粒子链（Alternating Gaussian Particle Chains）的**全区间相图、多粒子推广与高精度符号/区间算术验证库**。

## 背景与问题

考虑周期区间 $[0, N)$ 上的 $N=2k$ 个交替粒子（$A$ 类在偶数编号位，$B$ 类在奇数编号位），势函数为 $f(r)=e^{-r^2/4}$，周期化高斯核为

$$K_N(x) = \sum_{m\in\mathbb Z} e^{-(x + mN)^2 / 4}, \qquad E_N(x) = 1 - \frac{1}{N}\sum_{i,j=0}^{N-1} (-1)^{i+j} K_N(x_i - x_j).$$

在相邻粒子（含跨周期边界）满足最小间距约束 $\Delta x_i \ge h > 0$ 时，本项目完成了三项核心工作：

1. **$N=4$ 全区间 $h \in (0, 1]$ 相图与精确临界点**：
   - $0 < h \le h_{c} \approx 0.668019115487$：对称 `1+3` 构型 $(0, 2-h, 2, 2+h)$ 为全局极小；
   - $h = h_c \approx 0.668019115487$：发生一阶相变，三间隙锁死相 $(0, h, 2h, 3h)$ 能量低于对称构型（随后在 $h_{\text{instab}} \approx 0.668467328360$ 处对称构型失去局部 Hessian 稳定性）；
   - $h_c \le h \le 1$：三间隙锁死相 $(0, h, 2h, 3h)$ 成为全局极小，并在 $h=1$ 连续汇入等距构型 $(0, 1, 2, 3)$。
2. **推广到 $N=6, 8, 10$ 偶数粒子链**：
   - 发现小 $h$ 下不再是 `1+(N-1)` 聚类，而是倾向于分裂为两个紧致团簇（如 $N=6$ 的 `3+3`、$N=8$ 的 `3+5` 略优于 `4+4`）。
3. **高精度符号与区间算术严格验证（`mpmath` 55 位 + `mpmath.iv`）**：
   - 严格验证傅里叶稳定性界 $q = e^{-\pi^2/4} < 1/8$、尾项界 $\sum_{n\ge 2} n^2 q^{n^2} \le \frac{4q^4}{1-4q^2} \le \frac{q}{120}$、以及 $2q - 36\sum_{n\ge 2} n^2 q^{n^2} \ge \frac{17}{10}q$。
   - 补全了小间距局部证明中的包络定理消元步骤：须先由 $\partial_t E = 0$ 解出 $t^*(u,v) = u - v + O(\delta^3)$，再对 $u,v$ 求全导数得到 $\frac{dE}{du} = \frac{g''(\pi)-g''(0)}{2}v + O(\delta^3) > 0$。

## 快速运行

```bash
pip install -r requirements.txt
pytest -v
python scripts/run_all.py
```

运行 `scripts/run_all.py` 后将在 `figures/` 下生成：
- `figures/n4_phase_diagram.png`：$N=4$ 全区间能量曲线与间隙分岔相图
- `figures/multi_n_comparison.png`：$N=4, 6, 8, 10$ 最小能量对比图
