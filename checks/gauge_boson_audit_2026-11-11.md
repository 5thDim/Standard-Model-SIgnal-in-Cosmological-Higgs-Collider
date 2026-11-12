# Gauge boson signal 核查（2026-11-11）

## 结论与范围

按论文采用的约定，即常质量 massive vector、零化学势、实数质量指标
`nu = sqrt(m_A^2/H^2 - 1/4) > 0`，并只取 collapsed s channel 的 leading oscillatory nonlocal signal，
**当前 notebook 主流程与附录 B 的三个 C 一致，未发现这三个解析表达式的错误。**

核查对象为 `SMS in CHC.pdf` 第三章及附录 B、对应 TeX、`Gauge Boson Signal.nb` 和
`P2_reg0_analytic.wl`。同时核对本地 Qin–Xianyu 2023 原论文的式 (82)、(88)，与本文附录 B
式 (109)、(110) 的 double-folded scalar seeds 一致。

这次没有改动原始 TeX、PDF、notebook 或依赖脚本，也没有重新编译论文。
新增的 WLS 是复核脚本；在新的 Wolfram 15 kernel 中执行。它校验 notebook 的 SHA256 后选择
所需定义，跳过绘图、100 点扫描和旧归档中的绝对路径 `Get`，以免依赖已有 kernel 状态。

## 三个 C

令

\[
K=\frac{\Gamma^4(-i\nu)}{(4\pi)^{7/2}},\quad
I=\mathcal I_{\rm Bub}(-i\nu),\quad
P_1=\mathcal P^{(1)}_+(\nu),\quad
P_L=\mathcal P^{(2),L}_+(\nu),\quad
P_T=\mathcal P^{(2),+}_+(\nu).
\]

本次得到

\[
\begin{aligned}
C^{\rm Bub}&=K P_1^2 I,\\
C^{\rm Tri}&=\frac K3 P_1(P_L+2P_T)I,\\
C^{\rm Box}&=K\left[(P_L-P_T)^2\mathcal I^{LL}_{\rm Box}
 +\frac13P_T(2P_L+P_T)I\right].
\end{aligned}
\]

| PDF 公式 | 检查结果 |
| --- | --- |
| (103)，泡图 Gamma 闭式 | 与 `K P1^2 I` 的比值精确化简为 1；另检查 nu=0.01,0.1,0.6,1,2,5 |
| (107)，单个三角图 | 三个硬线极化的软环都精确等于 I/3，横向两 helicities 的因子 2 正确 |
| (125)，盒图 | 从 LL、LH、HH 的原始极化组求和得到相同表达式 |
| (126)–(127)，盒图五个角系数 | 五个多项式残差均精确为 0；公共 Gamma 归一化之比精确为 1；A3=-3 A2 成立 |

这里的 C 本身不包含图数、耦合和外部动量。第三章的 C^Tri 定义为一个 triangle placement，
两种 placement 的因子 2 放在总信号中；C^Box 定义为一个 representative box，collapsed s channel
两种 box 的因子 2 同样放在总信号中。与 notebook 比较时不能再把这些因子乘进 C。

## 独立软环计算

没有只比较优化前后的代码。对每条软线从空间极化张量重新开始：

\[
S_{ij}(\widehat q)=\delta_{ij}+(\beta-1)\widehat q_i\widehat q_j,
\qquad \beta=\frac{1+2x}{1-2x},\quad x=-i\nu.
\]

直接计算 `Tr[S(q) S(K)]`、`Tr[P_L S(K) S(q)]` 和
`Tr[P_L(k1) S(K) P_L(k3) S(q)]`，其中 `K=q+k_s`。
先用 Cartesian 单项式的方位角矩做积分，再把 z^r 用
`z=(K^2-q^2-k_s^2)/(2 q k_s)` 展开，独立归约到标量 Gamma 积分。
得到的 bubble、triangle-L 和 box-LL 相对于 notebook 的符号残差均为 0。

独立结果也给出

\[
\mathcal I^L_{\rm Tri}=\mathcal I^+_{\rm Tri}=\mathcal I^-_{\rm Tri}=I/3,
\quad
\mathcal I^{LH}_{\rm Box}=2I/3-2\mathcal I^{LL}_{\rm Box},
\quad
\mathcal I^{HH}_{\rm Box}=I/3+\mathcal I^{LL}_{\rm Box}.
\]

在三个质量指标及三个角度组合的交叉检查中，独立的六组软环结果相对 notebook 为数值零，
保留约 48 位有效精度。独立 z^r 展开与 J0–J4 的数值残差也为零。

## 时间积分、有限部和纵向接触项

- 从单顶点的 Gamma 积分重新计算 P1，与本文式 (34) 的比值精确为 1。
- 本文式 (115)–(117) 的 T_g、D_g 独立转录后，与 notebook 的 pMatrix 有两个精确零残差。
  regularized 3F2 的下参数导数包括 Gamma 分母的导数，不能用普通 3F2 的导数直接替代。
- `P2_reg0_analytic.wl` 的 inert-placeholder 方法保留 regulator 路径上的参数导数；
  实际执行后，有限部与优化主流程一致。
- 纵向源项分部积分也检查通过。对 SK 指标 a=±1，令
  `f(t)=t^(-3/2+i nu)(1+i a t)exp(-i a t)`，则
  `(d/dt+2/t)f=t^(-5/2+i nu)exp(-i a t)[(1/2+i nu)(1+i a t)+t^2]`。
  这给出式 (112) 中的三个 seed 组合及其正号。
- 不能删除式 (112) 的 `+P1/(nu^2+1/4)`。按本文 D 的定义，
  `d_tau1 d_tau2 D_s,++` 的等时分布项为 `+i H^2 tau^2 delta(tau1-tau2)`；
  从导数构造中扣掉它恢复 spatial longitudinal SK propagator，再乘 SK 因子，
  正好得到上面的正接触贡献。

在 nu=3/5，脚本直接给出

\[
\begin{aligned}
P_T&=0.496552665119014880418919820831
   +1.20229255009854174557228388109\,i,\\
P_L&=3.49769794788028052071618326871
   +6.83744965383156606331349582227\,i.
\end{aligned}
\]

该点 epsilon^(-3)、epsilon^(-2)、epsilon^(-1) 的系数数值均为零，约 73 位精度。
另从未先取有限部的原始 P2L、P2PM 出发，用 60 位输入、五个正 epsilon 值
`{1/1000,1/2000,1/4000,1/8000,1/16000}` 外推到零，得到相对解析 P 的最大误差：

| nu | 两个 P 中的最大相对误差 |
| --- | --- |
| 0.1 | 5.14e-20 |
| 0.6 | 4.00e-20 |
| 2 | 3.37e-19 |

## 发现的复现问题

1. **依赖脚本末尾的数值注释不适用于当前定义。** 第 65–66 行写的
   `2.870049... + 4.008603... I` 不是当前 `P2PMReg0[3/5]` 的结果，正确值见上。
   这属于注释错误；脚本实际有限部算法通过上述检查。
2. **旧归档的导入路径不可移植。** notebook 中的 `Get` 指向
   `D:\ChatGPT-Codex\Projects\Unprojected\2026-09-01\x\outputs\P2_reg0_analytic.wl`，
   且该 cell 目前为 `Evaluatable -> False`。因此当前优化主流程并不在运行时读取旁边的 WL 文件。
   若要重做原始有限部推导，需要显式加载当前项目中的 WL；本次已经这样做了。
3. **机器精度的小 epsilon 值会发生严重消去。** `P2PM[0.6,10^-8]` 含机器精度的 0.6，
   本次多次新 kernel 测试都明显偏离正确有限部；减小 epsilon 不能补救精度损失。
   应使用精确 `3/5` 或足够高精度的输入，再计算有限部或外推。
   `N[machine_result,50]` 不会恢复已损失的信息。
4. **旧归档的固定 epsilon 绘图只是近似。** 即使用高精度输入，epsilon=1e-5 仍有线性截断误差；
   例如 nu=0.6 时两个 P 的相对偏差约为几乘 1e-6。当前主流程使用精确有限部的高精度数值计算，
   不是这一固定 epsilon 近似。

## 复现与边界

在项目根目录执行：

```sh
/Applications/Wolfram.app/Contents/MacOS/WolframKernel -script checks/gauge_boson_audit_2026-11-11.wls
```

实际环境：Wolfram 15.0.1，Mac OS X ARM。数值测试使用精确质量指标或显式高精度输入。
本次核查结合了符号恒等式、独立软环推导、seed 文献对照及高精度 regulator 外推；
没有把振荡的原始四时间全积分另作直接数值积分。

本次读取版本的 SHA256：

```text
SMS in CHC.pdf
c5caed647fea732cb21371117475e9f9471e90fa288aa90adf07be740228de02
Gauge Boson Signal.nb
c232ddaf2a2ef7749f555458962db10832bfe04cdb7df2aa75270b0b313aa9e6
P2_reg0_analytic.wl
71b529994d5967cebfb334e72121cb743430232d522845d3606e8b051d5862d7
```
