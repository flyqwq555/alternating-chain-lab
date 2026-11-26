# $N=4$ 交替高斯粒子链相变的严格短证明

本文给出 $N=4$、周期 $L=4$、高斯势 $f(r)=e^{-r^2/4}$ 的最小间距问题中，**相变点 $h_c$ 的严格定位**与 **$t$ 消元（包络定理）步骤的完整推导**。

文中每一个数值常数都由 `scripts/verify_proof.py` 以 50 位十进制精度算出；所有级数截断误差用显式几何优级数控制（截断误差 $<10^{-30}$），非数值估计。符号 $\approx$ 后的数字为该脚本的输出。

---

## 0. 结论概览

| 编号 | 内容 | 状态 |
|---|---|---|
| 定理 1 | 两个候选构型的能量精确闭式 | 已证（初等） |
| 定理 2 | $F=E_A-E_B$ 在 $(0,1)$ 内唯一零点 $h_c$，且 $F'(h_c)>0$ | 已证（引理 6 的导数结构 + 机器认证） |
| 定理 3 | 1+3 拓扑的盒 $\delta\le u,v\le3\delta,\ |t|\le4\delta$ 内，$\delta\le0.0868$ 时极小构型唯一为 $(0,2-h,2,2+h)$ | 已证（含 $t$ 消元的完整推导） |
| 定理 4 | $h_c<h_{\text{instab}}$；两条候选能量在交点相等，下包络连续而导数跳变 | 已证（候选意义下的一阶相变） |
| 定理 5 | $0<h\le1/(2\pi)$ 时，全部合法构型的全局极小构型唯一为 $(0,2-h,2,2+h)$ | 已证（新增，见第 7 节） |
| — | $1/(2\pi)<h\le1$ 时**全部合法构型**下 $A$ / $B$ 的全局最优性 | **未证**，仅数值（第 6 节） |

---

## 1. 设定

周期化高斯核
$$K(x)=\sum_{m\in\mathbb Z}e^{-(x+4m)^2/4}.$$

**引理 1（傅里叶表示）** 令 $q=e^{-\pi^2/4}\approx 0.0848049724711137773$，$c=\sqrt\pi/2$，则
$$K(x)=c\left(1+2\sum_{n\ge 1}q^{n^2}\cos\frac{n\pi x}{2}\right).$$

*证明.* Poisson 求和：$f(r)=e^{-r^2/4}$ 的傅里叶变换为 $\hat f(\xi)=2\sqrt\pi e^{-\xi^2}$，周期 4 的周期化给出 $K(x)=\tfrac14\sum_n\hat f(\tfrac{\pi n}{2})e^{i\pi n x/2}$，整理即得。$\square$

**引理 2（对称性）** $K$ 为偶函数、4 周期，且关于 $x=2$ 对称，即 $K(2+y)=K(2-y)$。

*证明.* $K$ 偶且 4 周期，故 $K(2+y)=K(-(2+y))=K(y-2)=K(y-2+4)=K(2-y)$。$\square$

令 $\theta=\pi x/2$，$\delta=\pi h/2$，$g(\theta)=K(2\theta/\pi)=c\big(1+2\sum_{n\ge1}q^{n^2}\cos n\theta\big)$。
由引理 2，$g$ 偶、$2\pi$ 周期，且关于 $\pi$ 偶。于是
$$g'(0)=g'''(0)=0,\qquad g'(\pi)=g'''(\pi)=0. \tag{1.1}$$

四粒子能量（$s=(+1,-1,+1,-1)$）为
$$E=1-K(0)-\tfrac12\Big[K(x_0-x_2)+K(x_1-x_3)-K(x_0-x_1)-K(x_0-x_3)-K(x_2-x_1)-K(x_2-x_3)\Big]. \tag{1.2}$$

**引理 3** $g''(\pi)-g''(0)=4c\sum_{n\ \text{odd}}n^2q^{n^2}=\Delta\approx 0.300625807303069>0$。

*证明.* $g''(\theta)=-2c\sum n^2q^{n^2}\cos n\theta$，故 $g''(\pi)-g''(0)=-2c\sum n^2q^{n^2}\big[(-1)^n-1\big]=4c\sum_{n\ \text{odd}}n^2q^{n^2}>0$。$\square$

---

## 2. 定理 1：两个候选的精确能量

记
- **$A$（对称 1+3）**：$(0,\;2-h,\;2,\;2+h)$，周期间隙 $[\,2-h,\;h,\;h,\;2-h\,]$；
- **$B$（三间隙锁死）**：$(0,\;h,\;2h,\;3h)$，周期间隙 $[\,h,\;h,\;h,\;4-3h\,]$。

**定理 1** 对 $0<h\le1$，
$$E_A(h)=1-K(0)-\tfrac12\big[K(2)+K(2h)-2K(2-h)-2K(h)\big], \tag{2.1}$$
$$E_B(h)=1-K(0)-\tfrac12\big[2K(2h)-3K(h)-K(3h)\big]. \tag{2.2}$$

*证明.* 把 (1.2) 中的六个配对距离列出即可。对 $A$：同型对 $(0,2)$ 距离 $2$、$(1,3)$ 距离 $2h$；异型对距离为 $2-h,\;2+h\equiv 2-h,\;h,\;h$，其中用到周期 4 与 $K$ 的偶性。对 $B$：同型对均为 $2h$；异型对为 $h,h,h,3h$。代入 (1.2)。$\square$

---

## 3. 定理 2：相变点的严格定位

令
$$F(h):=E_A(h)-E_B(h)=\tfrac12\big[K(2h)+2K(2-h)-K(2)-K(h)-K(3h)\big]. \tag{3.1}$$

**引理 4** $F(0)=\tfrac12[K(2)-K(0)]<0$。

*证明.* 由引理 1，$K(0)-K(2)=2c\sum_n q^{n^2}\big(1-\cos n\pi\big)=4c\sum_{n\ \text{odd}}q^{n^2}>0$。$\square$

**引理 5** $F(1)=0$。

*证明.* $F(1)=\tfrac12[K(2)+2K(1)-K(2)-K(1)-K(3)]=\tfrac12[K(1)-K(3)]$，而 $K(3)=K(-1)=K(1)$（$K$ 偶且 4 周期），故为 $0$。$\square$

引理 5 是**退化现象**：$h=1$ 时两族构型重合为等距构型 $(0,1,2,3)$，此时可行集恰为单点，与相变无关。

**引理 6（导数结构）**
$$F'(h)=\tfrac12\big[2K'(2h)-2K'(2-h)-K'(h)-3K'(3h)\big], \tag{3.2}$$
$$F''(h)=\tfrac12\big[4K''(2h)+2K''(2-h)-K''(h)-9K''(3h)\big], \tag{3.3}$$
且 $F'(0)=F'(1)=0$。

*证明.* 前两式逐项求导。$F'(0)=\tfrac12[-2K'(2)]=0$（$K'$ 为奇函数）；$F'(1)=\tfrac12[-2K'(1)-K'(1)-3K'(3)]=\tfrac12[-3K'(1)+3K'(1)]=0$（因 $K'(3)=K'(-1)=-K'(1)$）。$\square$

**定理 2** 存在唯一 $h_c\in\big(0.6680191154868614433,\;0.6680191154868614435\big)$ 使 $F(h_c)=0$，且
$$F'(h_c)\approx 0.1009925510831>0.$$
从而
$$F(h)<0\ \ (0<h<h_c),\qquad F(h)>0\ \ (h_c<h<1).$$

*证明.* 分五段认证（每段的余量由引理 6 的导数界给出）：

1. $|F''|\le \tfrac12(4+2+1+9)\sup|K''|=8(\pi/2)^2\sup|g''|\approx 2.9742963$，
   $|F'''|\le 19(\pi/2)^3\sup|g'''|\approx 11.123037$，
   其中 $\sup|g^{(k)}|\le 2c\sum_n n^kq^{n^2}$（尾项由几何优级数控制，$<10^{-30}$）。
2. 在 $[0,10^{-3}]$ 上以步长 $2.5\times10^{-5}$ 认证 $F''>|F'''|\cdot s/2>0$，故 $F'$ 在 $(0,10^{-3}]$ 上为正（注意 $F'(0)=0$）。
3. 在 $[10^{-3},\,h^*-\rho]$ 上以步长认证 $F'>|F''|s/2>0$，其中 $\rho=2\times10^{-3}$。
4. 在 $[h^*-\rho,\,h^*+\rho]$ 内认证 $F''<0$；$F'$ 严格递减且两端异号，故 $F'$ 在该区间有唯一零点
   $h^*\approx 0.784532788955013$。
5. 在 $[h^*+\rho,\,1-10^{-3}]$ 上认证 $F'<0$；在 $[1-10^{-3},1]$ 上认证 $F''>0$ 且 $F'(1)=0$，故 $F'$ 在 $(1-10^{-3},1)$ 上仍为负。

于是 $F$ 在 $[0,h^*]$ 严格递增、在 $[h^*,1]$ 严格递减。由引理 4（$F(0)<0$）与引理 5（$F(1)=0$），$F$ 在 $(0,1)$ 内恰有一个零点 $h_c$，且 $h_c<h^*$ 处 $F'$ 为正（横截相交）。数值定位给出上面的区间。$\square$

**注** 定理 2 表明：若只比较 $A,B$ 两类构型，则 $h<h_c$ 时 $A$ 更低、$h\in(h_c,1)$ 时 $B$ 更低。

---

## 4. 定理 3：$t$ 消元与 1+3 拓扑内的唯一极小

固定粒子 0 于 $0$，另三个粒子记为
$$\pi+t-u,\qquad \pi+t,\qquad \pi+t+v,\qquad u,v\ge\delta>0.$$
由 (1.2)：
$$E(t,u,v)=1-g(0)-\tfrac12\big[g(\pi+t)+g(u+v)-g(\pi+t-u)-g(\pi+t+v)-g(u)-g(v)\big]. \tag{4.1}$$

### 4.1 方向导数

由 (4.1) 直接求导：
$$\partial_uE=-\tfrac12\big[g'(u+v)+g'(\pi+t-u)-g'(u)\big], \tag{4.2}$$
$$\partial_tE=-\tfrac12\big[g'(\pi+t)-g'(\pi+t-u)-g'(\pi+t+v)\big], \tag{4.3}$$
$$\partial_t^2E=-\tfrac12\big[g''(\pi+t)-g''(\pi+t-u)-g''(\pi+t+v)\big]. \tag{4.4}$$

**引理 7（$t$ 方向严格凸）** 在盒 $\mathcal B=\{\delta\le u,v\le3\delta,\ |t|\le4\delta\}$ 内
$$\partial_t^2E\ \ge\ \frac{a}{2}-\frac{M_4}{4}\big[t^2+(t-u)^2+(t+v)^2\big]
\ \ge\ cq\left(\frac{119}{120}-\frac{475}{8}\delta^2\right)>0,$$
其中 $a=g''(\pi)$、$M_4=\sup|g^{(4)}|$、$q=e^{-\pi^2/4}$、$c=\sqrt\pi/2$。

*证明.* 由 $g'''(\pi)=0$，对每个 $s$ 有
$$g''(\pi+s)=g''(\pi)+\int_0^s(s-r)g''''(\pi+r)\,dr,$$
故 $|g''(\pi+s)-g''(\pi)|\le M_4s^2/2$。代入 (4.4) 并合并公共项，得到
$$\partial_t^2E=\frac{a}{2}-\frac12\big[e(t)-e(t-u)-e(t+v)\big],$$
其中 $e(s)=g''(\pi+s)-a$，于是
$$\partial_t^2E\ge\frac{a}{2}-\frac{M_4}{4}\big[t^2+(t-u)^2+(t+v)^2\big].$$
在盒内 $t^2+(t-u)^2+(t+v)^2\le(16+49+49)\delta^2=114\delta^2$，并用
$$a=2c\sum_{n\ge1}n^2q^{n^2}(-1)^{n+1}\ge119cq/60,$$
$$M_4\le2c\sum_{n\ge1}n^4q^{n^2}\le25cq/12,$$
即得所示下界。这两个级数界与 (7.1) 相同。$\delta\le0.0868$ 时该下界仍大于 $0.54cq$。原先把点值 $\partial_t^2E(\delta,\delta,0)\approx0.074973$ 当成整个盒的下界并不成立。$\square$

### 4.2 $t$ 的存在唯一性与 $O(\delta^3)$ 精度

由 (1.1) 中 $g'''(\pi)=0$，对 $s$ 在 $0$ 附近展开
$$g'(\pi+s)=g''(\pi)\,s+\frac{g''''(\pi)}{6}s^3+R_\pi(s),\qquad |R_\pi(s)|\le\frac{\sup|g^{(5)}|}{24}|s|^4 . \tag{4.5}$$

**引理 8** 对每个 $(u,v)\in[\delta,3\delta]^2$，存在唯一 $t^*(u,v)\in[-4\delta,4\delta]$ 使 $\partial_tE=0$，且
$$\big|t^*(u,v)-(u-v)\big|\ \le\ 124.1\,\delta^3+306.8\,\delta^4 . \tag{4.6}$$

*证明.* 把 (4.5) 代入 (4.3) 并令其为零：
$$g''(\pi)\big(t-(u-v)\big)=\frac{g''''(\pi)}{6}\big[t^3-(t-u)^3-(t+v)^3\big]+R',\qquad |R'|\le3\cdot\frac{\sup|g^{(5)}|}{24}(7\delta)^4 .$$
在盒内 $\big|t^3-(t-u)^3-(t+v)^3\big|\le750\,\delta^3$，$\sup|g^{(5)}|\le0.1532466605$，$|g''''(\pi)|\le0.1488461012$，$g''(\pi)\ge0.1499461958$，即得 (4.6)。

严格凸性只给出驻点至多一个，存在性还需端点异号。记 $a=g''(\pi)$。在 $t=4\delta$ 处，对 $u,v\in[\delta,3\delta]$ 展开三项，
$$\partial_tE(4\delta,u,v)=-\frac12\big[a(u-v-4\delta)+R\big],\qquad |R|\le\frac{M_4}{6}\big[(4\delta)^3+(7\delta)^3+(7\delta)^3\big]=125M_4\delta^3.$$
因 $u-v-4\delta\le-2\delta$，故
$$\partial_tE(4\delta,u,v)\ge a\delta-\frac{125}{2}M_4\delta^3>0,$$
在 $\delta\le0.0868$ 时最后一项下界大于 $0.0068$。同法可得 $\partial_tE(-4\delta,u,v)<0$。结合引理 7，存在唯一 $t^*\in(-4\delta,4\delta)$。$\square$

### 4.3 包络定理：$u,v$ 方向的单调性（本文要点）

**定理 3（$t$ 消元）** 令 $\widetilde E(u,v):=E\big(t^*(u,v),u,v\big)$。则当
$$0<\delta\le\delta_0=0.0868\qquad\Big(\text{即 }0<h\le h_0=\tfrac{2\delta_0}{\pi}\approx0.055259\Big)$$
时，在 $\mathcal B$ 内严格有
$$\partial_u\widetilde E>0,\qquad \partial_v\widetilde E>0 .$$
从而 $\widetilde E$ 在 $\mathcal B$ 上的最小值唯一在 $u=v=\delta$ 取得。

*证明.* 因 $\partial_tE\big|_{t^*}=0$，有
$$\partial_u\widetilde E=\partial_uE\big|_{t^*}=-\tfrac12\big[g'(u+v)+g'(\pi+t^*-u)-g'(u)\big].$$
用 (4.5)（$g'''(0)=g'''(\pi)=0$）展开三项，并记 $t^*-u=-v+\eta$，其中 $|\eta|$ 由 (4.6) 控制：
$$\partial_uE\big|_{t^*}=\tfrac12\Delta\,v-\tfrac12 B,$$
$$B=g''(\pi)\eta+\frac{g''''(0)}{6}\big(3u^2v+3uv^2+v^3\big)+\frac{g''''(\pi)}{6}(t-u)^3+R .$$
在 $\mathcal B$ 内逐项估计（$\sup|g''''|\le0.1517797640$，$\sup|g^{(5)}|\le0.1532466605$）得
$$|B|\le b_3\delta^3+b_4\delta^4,\qquad b_3\approx31.89586067,\quad b_4\approx91.98630796 .$$
于是
$$\partial_u\widetilde E\ \ge\ \frac12\Delta\,v-\frac12\big(b_3\delta^3+b_4\delta^4\big)\ \ge\ \delta\left(\frac{\Delta}{2}-\frac{b_3}{2}\delta^2-\frac{b_4}{2}\delta^3\right),$$
其中 $\Delta/2\approx0.1503129037$。括号内为正等价于
$$\delta^2\le\delta_0^2,\qquad 0.1503129>15.9479\,\delta^2+45.9932\,\delta^3.$$
该粗下界的最大根约为 $0.0868207767454$，在根处下界等于零，不能推出严格正导数。固定取 $\delta_0=0.0868$ 时，括号内的余量约为 $6.87\times10^{-6}$，故严格不等式成立。由对称性（$u\leftrightarrow v$）同样得 $\partial_v\widetilde E>0$。$\square$

> **为什么必须先消元.** 若不对 $t$ 取驻值，(4.2) 的一阶形式是
> $$\partial_uE=-\tfrac12\big[g''(0)\,v+g''(\pi)\,(t-u)\big]+O(\delta^3),$$
> 其中 $g''(0)\approx-0.1506796$、$g''(\pi)\approx+0.1499462$。当 $t$ 未固定时该式**可取负值**：取 $\delta=10^{-4}$、$t=10\delta$、$u=5\delta$、$v=\delta$，得 $\partial_uE\approx-3.0\times10^{-5}<0$。
> 只有在 $t=t^*$（即 $t-u=-v+O(\delta^3)$）处，才有
> $$\partial_u\widetilde E=\tfrac12\big[g''(\pi)-g''(0)\big]v+O(\delta^3)=\tfrac12\Delta\,v+O(\delta^3)>0 .$$
> 这一步是原稿缺失的环节。

### 4.4 最终构型

由定理 3 取 $u=v=\delta$。此时由 (4.1) 与 (1.1)，
$$\partial_tE\big|_{(\delta,\delta,0)}=-\tfrac12\big[g'(\pi)-g'(\pi-\delta)-g'(\pi+\delta)\big]=0,$$
因为 $g'$ 关于 $\pi$ 为奇（$g'(\pi-\delta)=-g'(\pi+\delta)$）。结合引理 7 的严格凸性，$t=0$ 是 $t$ 方向唯一极小。

**推论** 当 $\delta\le\delta_0$ 时，在 1+3 拓扑且 $u,v\in[\delta,3\delta]$ 内，极小构型唯一为
$$(0,\ \pi-\delta,\ \pi,\ \pi+\delta)\quad\longleftrightarrow\quad(0,\ 2-h,\ 2,\ 2+h),$$
与定理 1 的 $E_A(h)$ 一致。特别地，原稿所述区间 $h<1/(500\pi)\approx0.00063662$ 远在认证范围 $h\le0.055259$ 之内。

---

## 5. 定理 4：候选意义下的一阶相变

**定理 4** 设 $h_{\text{instab}}$ 为对称构型 $A$ 失去横截（Hessian）稳定性的临界点，则
$$h_c\approx0.66801911548686144337\ <\ h_{\text{instab}}\approx0.66846732835982469912,$$
且 $A$ 在 $h_c$ 处的二阶变分为
$$\partial_s^2E_A\big|_{h_c}=-\tfrac12\big[K''(2)-2K''(2-h_c)\big]\approx2.2536\times10^{-4}>0 .$$

*证明.* $A$ 沿"大间隙分裂"方向 $s$ 的二阶变分为 $-\tfrac12[K''(2)-2K''(2-h)]$；其零点即 $h_{\text{instab}}$，由脚本区间定位。$h_c<h_{\text{instab}}$ 与二阶变分为正由同一脚本认证。$\square$

在交点处两条候选能量相等：
$$E_A(h_c)=E_B(h_c)\approx0.962490695544684293.$$
因此 $\min\{E_A,E_B\}$ 是连续函数，左右导数分别为
$$E_A'(h_c)\approx0.203226375782512,\qquad E_B'(h_c)\approx0.102233824699412.$$
这是候选意义下的**一阶相变**：能量连续，构型与能量导数发生跳变。$A$ 在 $h_c$ 的约束局部稳定性需要检查活跃约束的 KKT 条件，不能只用一个方向的二阶变分代替；现有计算也不能排除 $(h_c,h_{\text{instab}})$ 内存在其他稳定破缺分支。

---

## 6. 证明的边界（重要）

**已严格证明**：定理 1–3 的候选比较、$h_c$ 的唯一性与严格区间、$1+3$ 盒内的局部唯一极小，以及定理 4 中两条候选能量的交点与导数跳变。定理 5 进一步给出 $h\le1/(2\pi)$ 时全部合法构型的全局唯一性。

**未证明（仅为数值证据）**：

1. $1/(2\pi)<h<h_c$ 时 $A$ 是**全部合法构型**的全局极小；$h_c<h<1$ 时 $B$ 是全局极小。
   当前证据：全构型多起点 SLSQP 搜索（带严格可行性校验，见 `multi_particle.py`）在 $h\in(0,1]$ 的细网格上未发现更优构型。
2. $(h_c,h_{\text{instab}})$ 内不存在其他稳定对称性破缺分支。这需要额外的全局分支分析，不能由单个候选方向的曲率推出。

原稿中“$\mathcal D<5\delta^2$ 蕴含 $1+3$ 聚类”的定量缺口已在第 7 节补上。由此还得到一个小间距全局定理，但尚未覆盖到 $h_c$ 附近的全部区间。

---

## 7. 全局小间距定理

本节证明定理 5。全程使用角坐标
$$\theta=\frac{\pi x}{2},\qquad \delta=\frac{\pi h}{2},$$
于是 $0<h\le1/(2\pi)$ 等价于 $0<\delta\le1/4$。

### 7.1 核函数的一致估计

令
$$q=e^{-\pi^2/4}<\frac18,\qquad c=\frac{\sqrt\pi}{2},$$
$$T=\sum_{n\ge2}n^2q^{n^2},\qquad U=\sum_{n\ge2}n^4q^{n^2}.$$
由 $n^2\le4^{n-1}$、$n^4\le16^{n-1}$ 和几何级数，
$$T\le\frac{4q^4}{1-4q^2}\le\frac q{120},\qquad U\le\frac{16q^4}{1-16q^2}\le\frac q{24}. \tag{7.1}$$

对 $0<\theta<\pi$，利用 $|\sin n\theta|\le n\sin\theta$，
$$g'''(\theta)=2c\sum_{n\ge1}n^3q^{n^2}\sin n\theta
\ge2c(q-U)\sin\theta
\ge\frac{23}{12}cq\sin\theta>0. \tag{7.2}$$
因此 $g''$ 在 $[0,\pi]$ 上严格递增。

### 7.2 谱稳定性不等式

给定一个合法构型。令 $\Omega$ 为两段从正粒子沿正方向到下一个负粒子的圆弧之并，并记
$$J_n=\int_\Omega e^{in\theta}\,d\theta,\qquad
S_n=\sum_{j=0}^3(-1)^je^{in\theta_j}.$$
分部积分给出 $S_n=-inJ_n$，且由 (1.2)
$$E=1-\frac c2\sum_{n\ge1}q^{n^2}|S_n|^2.$$

旋转构型，使 $J_1$ 为非负实数。取半圆 $H=\{\cos\theta>0\}$，定义
$$\mathcal D=2-J_1=\int_{\Omega\triangle H}|\cos\theta|\,d\theta,
\qquad 0\le\mathcal D\le2. \tag{7.3}$$
令 $d=|\Omega\triangle H|$。由 $|\{|\cos\theta|<\varepsilon\}|\le2\pi\varepsilon$ 并取 $\varepsilon=d/(4\pi)$，得到
$$d^2\le8\pi\mathcal D. \tag{7.4}$$

对偶数 $n$，$J_n(H)=0$；对奇数 $n$，$J_n(H)=\pm2/n$。并且 $|\cos n\theta|\le n|\cos\theta|$。记 $b_n=J_n(\Omega)-J_n(H)$，则 $|b_n|\le d$，奇数 $n$ 时还有 $|\Re b_n|\le n\mathcal D$。于是
$$|S_n|^2-|S_n(H)|^2\le n^2(4\mathcal D+d^2)
\le(4+8\pi)n^2\mathcal D<36n^2\mathcal D. \tag{7.5}$$
第一模精确满足
$$|S_1|^2-4=-\mathcal D(4-\mathcal D).$$
记半圆能量
$$E_0=1-\frac{g(0)-g(\pi)}2.$$
将正负频率成对相加，并由 (7.1)，
$$\boxed{E-E_0\ge cq\mathcal D
\left(2-\frac{\mathcal D}{2}-\frac{18T}{q}\right)
\ge cq\mathcal D\left(\frac{37}{20}-\frac{\mathcal D}{2}\right).} \tag{7.6}$$
特别地，因 $\mathcal D\le2$，
$$E-E_0\ge\frac{17}{20}cq\mathcal D. \tag{7.7}$$

### 7.3 全局极小者的缺陷上界

合法间隙在固定平移后构成紧单纯形，故全局极小者存在。取试探构型
$$X_\delta=(0,\pi-\delta,\pi,\pi+\delta).$$
它的结构因子在奇数模为 $2\cos n\delta$，在偶数模为 $2(1-\cos n\delta)$，因此
$$\begin{aligned}
E(X_\delta)-E_0
={}&2c\sum_{n\ {\rm odd}}q^{n^2}\sin^2(n\delta)\\
&-2c\sum_{n\ {\rm even}}q^{n^2}(1-\cos n\delta)^2\\
\le{}&2c(q+T)\delta^2
\le\frac{121}{60}cq\delta^2. \tag{7.8}
\end{aligned}$$

设 $X$ 为任意全局极小者。由 $E(X)\le E(X_\delta)$、(7.7) 和 (7.8)，
$$\mathcal D\le\frac{121}{51}\delta^2<5\delta^2. \tag{7.9}$$
再代回较强的 (7.6)，当 $0<\delta\le1/4$ 时，
$$\mathcal D\le
\frac{121/60}{37/20-(5/2)\delta^2}\,\delta^2
\le\frac{968}{813}\delta^2
<\frac65\delta^2
\le\frac3{40}. \tag{7.10}$$

### 7.4 定量聚类引理

设某粒子到半圆 $H$ 最近端点的圆周距离为 $r$。若 $r\ge\delta$，该粒子左右各长 $\delta$ 的开弧内没有其他粒子，且 $H$ 在这两个开弧上取同一常值。由于 $\Omega$ 在粒子处跳变，其中一侧整段属于 $\Omega\triangle H$。距离函数为 $1$-Lipschitz，且 $\sin s\ge2s/\pi$ 对 $0\le s\le\pi/2$ 成立，所以
$$\mathcal D\ge\frac{2\delta}{\pi}(r-\delta).$$
若 $r<\delta$，该界自动成立。因此所有粒子到半圆端点的距离均不超过
$$R=\delta+\frac{\pi\mathcal D}{2\delta}. \tag{7.11}$$

**引理 9（$\mathcal D<5\delta^2$ 的聚类）** 若 $0<\delta\le0.0868$ 且 $\mathcal D<5\delta^2$，则四个粒子聚成两组，占据数为 $1+3$。

*证明.* 由 (7.11)，
$$R<\left(1+\frac{5\pi}{2}\right)\delta
\le\left(1+\frac{5\pi}{2}\right)0.0868<\frac\pi4. \tag{7.12}$$
两个中心长弧
$$[-\pi/4,\pi/4],\qquad[3\pi/4,5\pi/4]$$
均无粒子。在这两个长弧上 $|\cos\theta|\ge\sqrt2/2$。如果其中任一长弧上 $\Omega$ 与 $H$ 不一致，由于该长弧内没有粒子，$\Omega$ 在其上为常值，于是
$$\mathcal D\ge\frac{\pi}{2}\cdot\frac{\sqrt2}{2}
=\frac{\pi\sqrt2}{4}>\frac3{40},$$
与 $\mathcal D<5\delta^2\le5(0.0868)^2<3/40$ 矛盾。因此两个长弧上 $\Omega$ 与 $H$ 一致。

穿过每个粒子簇时，$\Omega$ 指示函数的净跳跃为 $+1$ 或 $-1$。粒子符号交替，所以每簇粒子数为奇数。四个粒子只能分成 $1+3$。$\square$

**注（仅有小 $\mathcal D$ 不足以进入 $3\delta$ 盒）** 取合法构型
$$(0,\ \pi-\delta,\ \pi+3\delta,\ \pi+4\delta),$$
其簇内参数为 $u=4\delta$、$v=\delta$、$t=3\delta$。第一模为
$$S_1=1+e^{-i\delta}-e^{3i\delta}+e^{4i\delta},$$
所以
$$\mathcal D=2-|S_1|
\le1-\cos\delta+\cos3\delta-\cos4\delta
\le4\delta^2<5\delta^2.$$
但它不满足 $u\le3\delta$。因此第 7.3 节的极小性强化 (7.10) 是必要的，不能只靠聚类引理把任意构型送入定理 3 的盒。

对定理 5 的全局极小者，用 (7.10) 代替 (7.12)：
$$R<\left(1+\frac{3\pi}{5}\right)\delta
\le\frac14+\frac{3\pi}{20}<\frac\pi4,
\qquad 0<\delta\le\frac14. \tag{7.13}$$
所以引理 9 的聚类结论在 $0<\delta\le1/4$ 仍成立。

### 7.5 聚类区域内的变分锁死

平移并循环重标号，把单粒子放在 $0$，另外三个写成
$$\pi+t-u,\qquad\pi+t,\qquad\pi+t+v,\qquad u,v\ge\delta.$$
由 (7.13)，设三粒子到 $\pi$ 的半径界为 $R$，则
$$|t|,\ |t-u|,\ |t+v|\le2R<\frac\pi2,
\qquad u+v\le2R<\frac\pi2. \tag{7.14}$$
两个外侧间距都大于 $\pi/2>\delta$，故全局极小者处 $t$ 可双向变化，满足 $E_t=0$。

直接求导得到两个精确积分恒等式
$$\begin{aligned}
E_u+E_t&=\frac12\int_0^v
\big[g''(\pi+t+s)-g''(u+s)\big]\,ds,\\
E_v-E_t&=\frac12\int_0^u
\big[g''(\pi+t-s)-g''(v+s)\big]\,ds.
\end{aligned} \tag{7.15}$$
在第一个积分中，$|t+s|<\pi/2$、$0<u+s<\pi/2$。由偶性、周期性和 $g''$ 在 $[0,\pi]$ 上严格递增，
$$g''(\pi+t+s)=g''(\pi-|t+s|)>g''(u+s).$$
第二个积分同理。因此 $E_u>0$、$E_v>0$。若 $u>\delta$，单独减小 $u$ 会降低能量且保持可行；若 $v>\delta$ 同理。故
$$\boxed{u=v=\delta}. \tag{7.16}$$

### 7.6 唯一驻点为 $t=0$

在 $u=v=\delta$ 时，傅里叶展开给出
$$E_t=c\sum_{n\ge1}nq^{n^2}(-1)^n(1-2\cos n\delta)\sin nt.$$
第一模为 $cq(2\cos\delta-1)\sin t$。其余模利用 $|\sin nt|\le n|\sin t|$ 控制。对 $0<|t|<\pi/2$，
$$\begin{aligned}
\frac{E_t}{\sin t}
&\ge cq\left(2\cos\delta-1-\frac{3T}{q}\right)\\
&\ge cq\left(\frac{39}{40}-\delta^2\right)
\ge\frac{73}{80}cq>0. \tag{7.17}
\end{aligned}$$
因此该区域内唯一的驻点是 $t=0$。结合 (7.16)，所有全局极小者均为
$$(0,\pi-\delta,\pi,\pi+\delta),$$
换回物理坐标即
$$\boxed{(0,\ 2-h,\ 2,\ 2+h)}. \qquad\square$$

**注** 定理 5 只覆盖 $h\le1/(2\pi)$。$h_c\approx0.668$ 两侧的全部合法构型分类仍需新的全局论证；第 6 节列出的未证部分没有被这里覆盖。

---

## 8. 复现

```bash
python scripts/verify_proof.py     # 逐项认证本文所有数值常数
pytest -v                          # 含 test_proof_constants.py
```

脚本使用的规则：

- 级数截断 $N=30$，尾项由 $\sum_{n>N}n^kq^{n^2}\le\sum_{n>N}n^kq^{(N+1)n}$ 控制，并用几何级数给出有限上界；
- 函数符号认证采用 Lipschitz 余量法：在网格点上验证 $f(x_i)>\sup|f'|\cdot s/2$，即保证整段成立；
- 工作精度 50 位十进制。

## 附：数值一览

| 量 | 值 |
|---|---|
| $q=e^{-\pi^2/4}$ | $0.0848049724711137773$ |
| $h_c$ | $0.66801911548686144337$ |
| $F'(h_c)$ | $0.1009925510831$ |
| $h^*$（$F'=0$） | $0.784532788955013$ |
| $h_{\text{instab}}$ | $0.66846732835982469912$ |
| $h_{\text{instab}}-h_c$ | $0.00044821287$ |
| $\Delta=g''(\pi)-g''(0)$ | $0.300625807303069$ |
| $\sup|g''''|$ | $0.1517797640$ |
| $\sup|g^{(5)}|$ | $0.1532466605$ |
| $b_3,\ b_4$ | $31.89586067,\ 91.98630796$ |
| $\delta_0$ / $h_0$（固定认证值） | $0.0868$ / $0.05525859624$ |
| 粗下界根估计 | $0.0868207767454$ |
| $1/(500\pi)$ | $0.0006366197724$（远在范围内） |
