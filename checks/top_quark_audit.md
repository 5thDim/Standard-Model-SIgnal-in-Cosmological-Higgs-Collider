# Top box：Aoki 两分量重写与复核

本轮将 `Top Quark Signal.tex` 的主线改为 Aoki 的两个独立 Majorana Weyl 场。已编译对应 17 页 PDF，未出现未定义引用、重复标签或 overfull/underfull 警告；构建记录见 `top_quark_build.txt`。物理条件仍为 exact de Sitter、恒定质量、零化学势及原来的轻 spectator Higgs 外腿。Gauge Boson 文稿及其偏振计算程序未修改。

## 与开始本轮时的结果比较

**最终领先非解析结果一致；本轮没有发现需要再改动的振幅系数。**

- 同分支最低阶 `ks^(3±4iμ)` 在四个 s-cut 外腿排列求和后消失。单图旋量迹一般不为零，单个硬块的完整 SK 积分也一般不为零。
- 首个振荡项为 `ks^(5±4iμ)`，保留正文五个角系数及 `A3=-3 A2`。硬传播子动量导数与两个软端点的 first descendants 缺一不可。
- 异分支给出各向同性的非振荡 `ks^3`。它对软动量分量非解析，不能作为多项式接触项删除。
- 正文给出完整 `Nc gt^4 H^4`、外腿分母、硬尺度、软 Gamma 卷积及 SK 硬 seed 的归一化；`Nc=3` 只计一次。
- 固定非零 μ、中心化硬动量下，下一阶分别为 `ks^(7±4iμ)` 与 `ks^5`。这是逐分支的渐近精度；不声称在质量趋零或角系数为零附近有统一的相对误差界。

旧检查记录提到的 regularized `3F2` 分母 Gamma 因子遗漏，在本轮开始前的工作区版本中**已经修复**。本轮核对了 Qin–Xianyu 原 PDF 的 (92)、(93)、(138)，并保留修正；不能把这件历史修复列为本轮新发现。

## 本轮补齐和修正的实质内容

1. **真正采用 Aoki 变量，而非只把 Dirac 分量称为 Weyl。** 从 `t=(ψ1+iψ2)/sqrt(2)`、`tc=(ψ1-iψ2)/sqrt(2)` 推得两份含 `1/2` 的 Majorana 作用量。说明 `tc` 的电荷共轭意义，以及忽略规范相互作用时逐颜色分解的适用性。
2. **全部正常/反常收缩与顶点。** 固定 `epsilon_lower=i sigma2`、`epsilon_upper=-epsilon_lower`、所有场统一的 Fourier 符号和有符号 lesser 收缩。反向线同时交换端点、类型、转置自旋指标、反转动量并带 Grassmann 负号。共轭与反向是不同操作。给出两种 Yukawa 顶点，既可保留双线性 `1/2` 枚举槽位，也可使用已消去 `2!` 的顶点，不能同时用两种计数。
3. **独立 Majorana 计数。** 八场 Wick 定理的 105 个配对中有 48 个连通单圈配对，即 3 个无向周期各 16 种槽位连接。每顶点 `1/2` 后，每个无向周期权重为 1；写成 6 个定向周期时每个为 `1/2`。两份 Majorana 恢复 Dirac 结果，再乘颜色数。没有为了匹配旧振幅调整因子。
4. **mode 字典。** Aoki 的 Whittaker `uλ,vλ` 与旧 Hankel `F,G` 的精确共同相位为 `exp(3πi/4)`；`vλ` 不是旧四分量负频模式。全文恢复 H，使用 `μ=mt/H`。
5. **首个存活项的完整阶数检查。** 加入精确 `0F1` 软端点表示；同时保留硬动量、外腿能量移动及软 descendants。外腿交换抵消、SK 晚时端点抵消和单图迹是三种不同事项。
6. **修正验证强度的表述。** 两种 scalar-seed 表达式互相比对属于一致性检查；有限时间区间的展开验证不等于完整 BD 振幅验证。新增实际物理质量指数下的直接 BD 时间积分，而非只检验 conformal-scalar 特例。

## 参考文献与组织

- [Aoki et al., arXiv:2605.28054](https://arxiv.org/abs/2605.28054)：本地 `ref/` PDF，重点 (2.5)–(2.29)。仅使用自由场和两分量约定，未将其 bubble 振幅或 Yukawa 三点函数 cancellation 当成 box 结果。
- [Qin–Xianyu, arXiv:2301.07047](https://arxiv.org/abs/2301.07047)：scalar seed 及 regularized 超几何函数。
- [Qin–Xianyu, arXiv:2304.13295](https://arxiv.org/abs/2304.13295)：双软区域因子化。
- Gauge Boson 文件提供“硬子图—旋量/偏振缩并—标量卷积—角系数”的组织方式；其零阶硬线替换不足以计算这里的首个振荡项。

正文现在按 setup → 两分量传播子和规则 → Majorana box Wick 缩并 → collapsed 展开和抵消 → scalar 积分 → 最终信号组织。仅在附录保留四分量 Dirac 对照。正文的 `R^{rs}` 是由正常/反常两分量收缩定义的质量投影，每个元素仍为 `2×2`，不代表额外物种。

## 新增程序和验证层级

| 文件 | 检验内容 | 独立性与范围 |
|---|---|---|
| `aoki_majorana.py` | Whittaker 两分量收缩、Wick 枚举、Dirac 对照 | 两分量从 mode operator 展开构造；Dirac 独立用 Hankel |
| `aoki_majorana_checks.py` | EOM、BD、归一化、全部 SK 和共轭、精确软分支、48 配对、完整圈收缩 | 覆盖 μ=0.2,0.73,2,5；圈图覆盖前三个质量和全部 16 种 SK；generic edge momenta 的局部代数检查 |
| `aoki_collapsed_checks.py` | 四个 s-cut 排列、完整软端点、硬动量及外腿能量移动 | μ=0.2,0.73,2，两种软动量几何，三个软尺度；有限时间区间积分，非 BD 振幅 |
| `aoki_time_integral_checks.py` | 直接 Hankel 时间积分及动量导数对照 scalar seed | same-SK ordered triangles 与 opposite-SK 分别按 BD 旋转；物理 μ=0.2,0.73，含 mixed；共同 δ=0.8；同时验证 P1/P2，不是 δ→0 极限 |
| `top_box_signal_checks.py` | 最终公式的两端 Bose 对称、pair 交换、宇称及量纲 | 组装公式的一致性检查，不是独立圈积分 |
| `top_box_signal.py` | 从守恒的四个三维动量计算最终 clock、mixed 和总和 | 包括两份 Majorana 和 Nc；返回软展开参数；只是结果求值器，不是独立验证 |

保留并重新运行：

- `weyl_signal_checks.py`：模式/传播子局部恒等式、Pauli trace、Gamma 卷积、五角系数及 mixed 因子。
- `weyl_hard_expansion_check.py`：原 Hankel 表示的有限区间硬展开检查。
- `weyl_seed_checks.py`：初等 conformal-scalar 积分、regularized seed、共同 regulator 路径、精度和 mixed 步长稳定性。
- `weyl_scalar_seed.py`：经验证后复用的 scalar seed，不包含 Aoki bubble 振幅。已更新文档说明验证层级。

每个检查脚本对应的 `.txt` 保存实测结果。代表性结果：

- Whittaker/Hankel、EOM 和归一化残差小于 `3e-38`；正常/反常 SK 与 Dirac 映射小于 `8e-16`。
- 显式八槽位 Wick 求和与两份 Majorana/Dirac 完整圈图的一致性残差小于 `4e-14`（双精度矩阵求和）。
- 精确 `0F1` 软端点与 Hankel 的 Bessel 分支拆分残差小于 `3e-41`。
- 双软展开的相对残差在软尺度减半时约缩小四倍；最小尺度的误差为 `7e-6` 到 `2.2e-5`。有限区间 N=32→48 的系数变化约 `2.8e-3`，与带时间排序 cusp 的张量求积一致；不据此声称无限时间积分精度。
- 直接 BD 时间求积的三个配置：N=48→72、积分尾端26→30，最终对照 seed 的相对差分别约 `2.7e-9`、`6.7e-10`、`1.9e-6`；求积变化分别约 `1.4e-8`、`3.5e-9`、`1.3e-5`。clock 的 H11 使用首个非零系数所需的额外时间幂。直接 regulated P1/P2 在 μ=0.2,0.73 的误差分别约 `2.6e-9`、`1.0e-9`，其动量导数由 mode 方程计算，不调用 scalar bootstrap。
- hard 系数在 22/35 位精度间的变化小于 `1.2e-21`；共同 regulator 三条路径外推残差小于 `1.1e-8`。
- μ=0.73 的 mixed seed 步长减半变化约 `1.9e-15`。

μ=0.73 时仍得到：

```text
P1 = -1.59813070834661283 - 2.34621969383716141 i
P2 =  1.92388500855456824 + 1.95859683065576684 i
H11(mixed) = -1.3885002801988139277
```

## 复现

在项目根目录执行（需要 numpy、mpmath、sympy 和 TeX Live/latexmk）：

```bash
python3 checks/aoki_majorana_checks.py
python3 checks/aoki_collapsed_checks.py
python3 checks/aoki_time_integral_checks.py
python3 checks/weyl_signal_checks.py
python3 checks/weyl_hard_expansion_check.py
python3 checks/weyl_seed_checks.py
python3 checks/top_box_signal_checks.py
python3 checks/top_box_signal.py --mu .73 --soft .02
latexmk -pdf -interaction=nonstopmode -halt-on-error 'Top Quark Signal.tex'
```

最后一个 Python 命令使用 H=gt=1、Nc=3、kL=(0.8,0,0.6)、kR=(0.3,1.2,-0.4)、ks=(0,0,0.02) 的示例；结果见 `top_box_signal_example.json`。这些是求值器参数，不是另加的物理假设。

## 尚未独立完成的检查

- 未对未展开的完整 box 同时做四个时间及三维 loop momentum 的数值积分。
- 直接 BD-contour quadrature 在共同正 regulator 下验证了 scalar seed；其 δ→0 的最终振幅仍由解析 seed 求和、不同 regulator 路径和精度检查确认，没有第二套直接无 regulator 数值求值。
- 未用 Aoki 的完整 spectral/Mellin–Barnes 方法重新计算四顶点 box；其 bubble 结果不提供这项验证。
- 解析接触项的重整化、重叠 collapsed 通道和转为曲率 trispectrum 不属于此处已交付的领先非解析项。
