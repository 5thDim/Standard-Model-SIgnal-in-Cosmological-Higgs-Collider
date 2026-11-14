# Signal Plots

打开根目录的 `Signal Plots.nb` 即可查看原有十张图、六张曲率 tNL 角分布预览、两张质量依赖图、两张 trispectrum 对比图及一张扩展范围的 T/E 图。正常运行只加载本目录中的六个模块，不读取原 notebook、原辅助脚本或论文。原文件均保留。

| 文件 | 职责 |
|---|---|
| `coefficients.wl` | 最终解析有限部、特殊函数稳定求值、复数系数接口及总角向组合 |
| `observables.wl` | 第五节的 W/Z/top 复系数组装、曲率 clock amplitude 与相位 |
| `plot_data.wl` | 精确默认网格、复数扫描、固定质量扫角、幅度/相位投影、WXF 保存与缓存 |
| `plot_style.wl` | Latin Modern Roman 字体、字号、线宽、图尺寸、留白、图例、渐变色及分段色条 |
| `nu_scan.wl` | 根据已算好的复数扫描画系数图及曲率 tNL/相位质量图；相位跳变处断线 |
| `angular_maps.wl` | 根据已算好的角网格画幅度和相位填色等高线地图 |
| `build_notebook.wl` | 在新内核中生成并求值简短 notebook，保存单图及 2×2 组合预览 |
| `validate.wl` | 可选的只读移植回归；读取原 notebook 的解析定义及第 8 节实例作比较 |

## 修改参数与重画

参数集中在新 notebook 第 1 节的第二个输入格：

- `gaugeNuGrid`：`Range[1/100,1,1/100]`，100 点。
- `topNuGrid`：0.1–1 步长 0.01，1.05–5 步长 0.05，共 171 点。
- `gaugePrecision=40`、`topPrecision=50`、`angularPrecision=60`。
- `angularNu=1`、`phi3=0`、`cosineGrid=Subdivide[-1,1,200]`，即 201×201 网格。
- `amplitudeUnit=10^-6`：色条标为“|Ct| × 10^6”，刻度 5 对应实际幅度 5×10^-6。幅度未取对数或归一化。
- `magnitudeEdges`、`phaseEdges` 为原第 8 节色阶。更换二维参数后可改为 `Automatic`。

从上到下运行会自动使用匹配的 `data/*.wxf`；参数、精度或 `coefficients.wl` 的哈希改变时重算。`recomputeData=True` 可强制重算。只改样式时，修改并重新 `Get` 加载 `plot_style.wl`，再运行第 3–5 节即可，不用再算特殊函数。

Gauge 的七条曲线依次使用紫、蓝、青、绿、深金黄、橙、红；在 `plot_style.wl` 的 `GaugeScanColors` 中调整。两张 Gauge 图使用相同配色。

Top 的五条曲线依次使用紫、蓝、绿、深金黄、红；在 `plot_style.wl` 的 `ScanColors` 中调整。两张 Top 图使用相同配色。

第 5 节另含 ϕ3=π/2、π 的四张图：2×2 排列，上排 π/2、下排 π，左列幅度、右列相位；仍为 ν=1、201×201 网格。对应数据分别保存在 `top_angular_phi_pi_over_2.wxf` 与 `top_angular_phi_pi.wxf`。这一组的色阶在该小节的 `additionalMagnitudeEdges`、`additionalPhaseEdges` 中设置，并在两行间统一。

所有 notebook 路径从 `NotebookDirectory[]` 构造；两个命令行脚本从自身所在目录定位工程。没有绝对路径或 Windows 路径。命令行可在项目根目录运行：

```sh
wolframscript -file plots/build_notebook.wl
wolframscript -file plots/validate.wl
```

## 解析式与接口

加载模块后：

```wl
g = GaugeCoefficients[3/5, 60];  (* CBub, CTri, A1, A2, A3, B, D *)
t = TopCoefficients[1, 60];     (* A1, A2, A3, B, D *)
c = AngularCoefficient[t, 1/3, -1/2, 0];
{Abs[c], Arg[c]}
```

返回值均为字符串键的 Association；`"D"` 不覆盖内置 `D`。`GaugeCoefficientExpressions[nu]`、`TopCoefficientExpressions[nu]` 返回对应的精确解析表达式。请优先使用分数等精确输入；机器数如 `0.6` 会被拒绝，避免先丢失精度再增加精度。

Gauge 中 `nu=Sqrt[(mA/H)^2-1/4]>0`；top 中 `nu=mt/H>0`，即论文的 fermion 质量指标。沿用原正频软支、相位、复数平方以及单个代表拓扑的归一化；没有加入额外耦合、图形重数或动量因子。

Gauge 的有限时间积分移自最终显式系数节的 `tFiniteCore`、`dFiniteCore`，保留纵向接触项；组装方式等同于主流程 `commonLoop`、`allTemplates`。Top 使用最终 `topPReduced` 的有限表达式及主流程 `softBTop`、`coefficientTemplate`。没有固定 epsilon 近似或重新推导调节参数极限。

Gauge 保留四个超几何母函数，包括 regularized 3F2 的参数导数。混合参数导数使用原 notebook 的 Cauchy 轮廓求值，采样数随精度增长；这是解析导数的数值计算，轮廓半径不是调节参数 epsilon。其余特殊函数继续由 Wolfram 的原解析延拓求值。

总角向函数始终先计算

```wl
A1 + A2 (x^2+y^2) + A3 x^2 y^2
 + 4 B x y Sqrt[1-x^2] Sqrt[1-y^2] Cos[phi3]
 + D (1-x^2) (1-y^2) Cos[2 phi3]
```

再求 `Abs` 或 `Arg`。Gauge 的角向函数仅代表 box，`CBub` 和 `CTri` 保留各自定义。

## 曲率扰动的 clock observables

`observables.wl` 使用论文第五节的归一化，并保留 `coefficients.wl` 中的原始定义。
默认参数为 Lu et al. (2019), [arXiv:1907.07390](https://arxiv.org/abs/1907.07390)
第 2 节和图 4 的取值，其中 Higgs 背景按本项目要求取 `h0/H=1`：

```wl
clockParameters = Lu2019Parameters[];
(* <|"g"->3/5, "cW"->22/25, "yt"->1, "lambda"->1/100,
     "Rh"->7/50, "Pzeta"->2/10^9, "h0OverH"->1|> *)

zClock = ClockCoefficients["Z", 1, 60, clockParameters];
wClock = AssembleClockCoefficients["W", zClock["RawCoefficients"], clockParameters];
topClock = ClockCoefficients["t", 1, 60, clockParameters];
zValues = ClockObservables[zClock, Pi/2, Pi/2, 0];
(* keys: Species, C, ComplexAmplitude, tNL, Phase *)
```

Notebook 第 1 节的新增小节集中设置 `clockParameters`，并给出三种粒子的参考数值表。
示例沿用 `angularNu`，取 `theta1=theta3=Pi/2` 及现有 `phi3`；该点用于展示调用，尚未进行角度寻优。
`nu` 是各粒子自己的质量指标，W/Z 对应 `m/H=Sqrt[nu^2+1/4]`，top 对应 `m/H=nu`。
`lambda` 用于记录参考参数，本步的转换公式不直接依赖它，也不据此自动固定粒子质量。

`AssembleClockCoefficients` 接受已算好的原始系数，生成两个五项角向 Association：

- `AngularCoefficients` 表示第五节的 `C_W`、`C_Z` 或 top 的 `C_t`。
- `ComplexAmplitudeCoefficients` 包含曲率归一化和整体符号，可直接用于后续角网格。

W 的组合是 `CBub-2 g^2 (h0/H)^2 CTri+g^4 (h0/H)^4 CBox`；Z 在组合中取
`g->g/cW`。W 的振幅包含两个电荷态，top 使用 `3 yt^4 Rh^4/(32 Pi^2 Pzeta)`，
相位为 `Arg[-Ct]`。`ClockObservables` 的角度输入均为弧度，返回的 `C` 保留复数，
并在复数相加之后计算幅度和相位。
`ComplexAmplitude=tNL Exp[I Phase]`，因此余弦形式等于 `Re[ComplexAmplitude Exp[I oscillationPhase]]`；
这里不再额外乘 2。振幅严格为零时，相位返回 `Missing["ZeroAmplitude"]`。

同一质量指标下，W 和 Z 可共享原始 gauge 系数；改变外动量方向时复用组装结果，无需重算特殊函数。
修改基准参数后应重新组装；原始复系数缓存仍可复用。原有系数图和 WXF 数据不因本步而改变。

新增组装层已在 Wolfram 15.0.1 中通过 54 项检查，包括复数拓扑求和、W/Z 重数与耦合权重、
top 相位负号、原四点函数到曲率 trispectrum 的归一化及零振幅处理。
Notebook 的三个新增输入已在新内核中重放，60 位精度的参考数值与保存的表格一致。

## 六张 tNL 角分布预览

Notebook 第 6 节绘制 Z 和 top 的实际无量纲 `tNL`，沿用第 1 节的 `angularNu=1`、
`angularPrecision=60`、`clockParameters` 和 `clockAssemblies`。
横纵坐标为 `theta1`、`theta3`，均从 0 到 Pi，默认各 201 点。
三列分别为 `phi3=0,Pi/2,Pi`；Z 与 top 分别使用自己的统一色阶，使同一粒子的三张图可直接比较。
色条使用实际 tNL 数值，未额外缩放、取对数或按最大值归一化。
色条刻度采用科学计数法，标题使用上下对齐的 `t_NL^Z`、`t_NL^t`。
粒子种类、nu 和 phi3 放在图内右上角的白底黑色直角边框中。
PNG 使用不透明白底和黑色文字，避免显示主题影响坐标及色条的可读性。

预览位于 `previews/`：

| 粒子 | phi3=0 | phi3=Pi/2 | phi3=Pi |
|---|---|---|---|
| Z | `z_tnl_phi_0.png` | `z_tnl_phi_pi_over_2.png` | `z_tnl_phi_pi.png` |
| top | `t_tnl_phi_0.png` | `t_tnl_phi_pi_over_2.png` | `t_tnl_phi_pi.png` |

`tnl_angular_2x3.png` 为两行三列总览，上排 Z、下排 top。
`tnl_angular_summary.txt` 记录各图的实际幅度范围和数值核验。
本阶段只更新预览，不导出插图 PDF，也不改动论文。

修改第 6 节的 `clockPhiCases`、`clockThetaGrid`、`clockBandCount` 可以调整角度切片、网格和色阶数。
改变 nu 或物理参数后，先重新运行第 1 节的组装，再运行第 6 节。
仅改配色、字号等样式时，重新加载 `plot_style.wl` 和 `angular_maps.wl`，
再运行第 6 节的绘图和导出两格即可。

新增接口 `GenerateClockAngularData` / `CachedClockAngularData` 位于 `plot_data.wl`，
`ClockAngularMapPlot` 位于 `angular_maps.wl`，坐标与标签样式由 `plot_style.wl` 的
`ClockMapPlotOptions` / `ClockMapLegendLabel` / `ClockMapColorLegend` 提供。
新缓存保留归一化后的复数振幅，并检查物理参数、nu、phi3、theta 网格、精度及组装/源码签名；
文件名为 `data/{z,t}_tnl_angular_phi_*.wxf`，与原始系数角图分开。

## 曲率 tNL 与相位的质量依赖

Notebook 第 7 节固定 `clockMassAngles={Pi/2,Pi/2,0}`，顺序为 `theta1,theta3,phi3`。
沿用角分布图所选的构型，并在质量扫描中保持这一构型不变。
采用第 1 节的 `clockParameters`，保持 `h0/H=1`，独立扫描各粒子的质量。
横轴为各粒子自身的 `m/H`，质量输入使用精确数：

- W、Z：`0.51..1` 步长 `0.01`，`1.025..3` 步长 `0.025`，共 130 点。
- top：`0.1..1` 步长 `0.01`，`1.025..3` 步长 `0.025`，共 171 点。
- `clockMassPrecision=60`。W/Z 使用 `nu=Sqrt[(m/H)^2-1/4]`，top 使用 `nu=m/H`。

W/Z 仅画 `m/H>1/2` 的振荡分支，未取阈值点。
三条曲线先相加复系数，再取幅度或相位；W、Z、top 在两图中分别使用红、蓝、黑色。
三条曲线均为实线。
两张质量图不显示固定角度标记；三种粒子的红、蓝、黑线图例统一放在图内的白底黑色直角方框中。
Notebook 中保存压缩 PNG 预览，重新运行绘图格可生成原生图形。
`previews/clock_tnl_vs_mass.png` 显示实际无量纲 tNL，纵轴为对数刻度。
`previews/clock_phase_vs_mass.png` 显示 `(-Pi,Pi]` 内的主值相位，单位为弧度，
在相位跨越支切处断线，不展开相位。top 沿用 `Arg[-Ct]`。

原始数据存为 `data/gauge_clock_mass_scan.wxf` 与 `data/top_clock_mass_scan.wxf`；
W/Z 共享同一份 gauge 扫描，各自按耦合权重和重数组装。
`AssembleClockMassScan` 位于 `observables.wl`，不重新调用特殊函数。
`data/clock_mass_observables.wxf` 保留三种粒子的质量、nu、C 和归一化复振幅，
并记录精度、参数、固定角度和原始系数源码签名。
`ClockMassPlotData` 位于 `plot_data.wl`；`ClockMassScanPlot` 位于 `nu_scan.wl`，
只从保存的复数数据提取幅度或相位。

改变固定角度或物理参数后，运行第 7 节的数据组装及绘图、导出格即可复用原始系数缓存。
改变质量网格、精度或系数源码会重新计算相应的独立质量扫描缓存。
只改样式时重新加载模块后运行两个绘图格及导出格。
本阶段只输出 PNG 预览，不改动论文及插图 PDF。

## Trispectrum 与去共同包络的对比

Notebook 第 8 节绘制带正负号的振荡非局域 collapsed 模板 `T_NS`，以及 `T_NS/E`。
动量比为 `r1=ks/k12`、`r2=ks/k34`，其中 `k12=k1+k2`、`k34=k3+k4`。
横纵坐标为 `x=-log10(r1)`、`y=-log10(r2)`，默认网格为 `Subdivide[0,2,200]`，共 201×201 点。
固定 `theta1=theta3=Pi/2, phi3=0`，其余参数沿用第 1 节，`lambda=1/100`。

`Lu2019MassBenchmark` 按 Lu (2019) 的式 (11)、(45) 使用环修正主导的 W/Z 质量，
top 则使用 `mt=yt h0/Sqrt[2]`。三种质量约为 `(0.557624,0.633664,0.707107) H`。
`AssembleClockShape` 保存该构型的复振幅、质量指标、相位和原始组装系数；
`ClockCollapsedSignal` 用它们计算单个动量比构型，不再求值特殊函数。

使用的关系为：

```text
E = (r1 r2/(r1+r2))^3
Psi_I = 2 nu_I Log[r1 r2/4] + phase_I
T_NS/E = tNL_W Cos[Psi_W] + tNL_Z Cos[Psi_Z] + tNL_t r1 r2 Cos[Psi_t]
T_NS = E (T_NS/E)
```

两张图均保留信号符号，采用以零为中心的线性色条、科学计数法刻度和相同的红蓝配色；
各图的色条范围独立设置，数值未取绝对值或按最大值归一化。
`T_NS/E` 只去掉共同包络，top 额外的 `r1 r2` 衰减仍然保留。
`r_i` 接近 1 的区域属于 collapsed 模板的外推；深 soft 区域才用于解释 collapsed 近似下的信号。

- `previews/clock_trispectrum.png`：实际 `T_NS`。
- `previews/clock_trispectrum_over_envelope.png`：实际 `T_NS/E`。
- `previews/clock_trispectrum_summary.txt`：物理参数、数值范围及核验结果。
- `data/clock_shape_log_grid.wxf`：60 位精度的 `Envelope`、`T`、`TOverE` 和三个 `SpeciesReduced` 矩阵；
  行对应 y，列对应 x。文件同时保留精确坐标网格、参数、构型、复系数及源码签名。

`GenerateClockShapeData` / `CachedClockShapeData` 位于 `plot_data.wl`；
`ClockShapeMapPlot` 位于 `angular_maps.wl`。修改 `clockShapeLogGrid` 可调整范围和分辨率，
修改 `clockShapeBandCount` 可调整色阶数（取正偶数，至少 2）。
只改样式时重新加载样式和地图模块，再运行第 8 节的两个绘图格及导出格即可。
Notebook 保存压缩 PNG 预览；重新运行绘图格会生成原生图形。
本阶段只更新 notebook 和 PNG 预览，不改动论文及插图 PDF。

## 扩展到 0..5 的 T/E 图

Notebook 第 9 节单独绘制 `T_NS/E`，两轴均取 `-log10(r_i)=0..5`，即 `r_i=10^-5..1`。
保持 201×201 网格、60 位数据精度和第 8 节的质量、参数、构型及样式。
本图使用数学字体的花体 T 字形，并通过 `"ShowParameters"->False` 隐去 lambda、h0 和角度标记。
复用第 8 节的 `clockCollapsedShape`，或从原 `clock_shape_log_grid.wxf` 中加载其 `Shape`。
本步不重新求值特殊函数，也不计算或绘制新的 T、E 网格。

固定构型下 `T_NS/E` 仅依赖于 `x+y`，`GenerateReducedClockShapeData` 因而只求值
401 个不同的坐标和，再填入二维网格。`CachedReducedClockShapeData` 单独缓存此结果。
top 的额外 `r1 r2` 因子仍然保留；物种分量仅以一维 `SpeciesReducedSeries` 保存。
两个接口位于 `plot_data.wl`，绘图仍使用 `ClockShapeMapPlot`。

- 预览：`previews/clock_trispectrum_over_envelope_0_to_5.png`。
- 数据：`data/clock_shape_reduced_log_grid_0_to_5.wxf`。
- 核验：`previews/clock_trispectrum_over_envelope_0_to_5_summary.txt`。

修改 `clockWideLogGrid` 可调整范围及网格（须均匀）；只改样式时重画第 9 节即可。
原 0..2 的两张图、数据及 notebook 小节均保留，论文和插图 PDF 不变。

## 数据结构

WXF 保留任意精度复数，绘图时只将显示副本转换为机器精度。

| 数据 | 主要字段 |
|---|---|
| ν 扫描 | `Model`、`Precision`、精确 `NuGrid`、`CoefficientNames`；`Rows` 每项为 `<|"Nu"->nu,"Coefficients"-><|系数名->复数值|>|>` |
| 二维扫角 | `Model`、`Nu`、`Phi3`、`Precision`、精确 `CosGrid`、`GridDimensions`；`Coefficients` 保存仅计算一次的系数；`Rows` 为 `{x,y,complex C}`，y 外层、x 内层 |

两者还保存 `SchemaVersion` 和 `CoefficientSourceHash`。`NuPlotData` 与 `AngularPlotData` 只转换已有复数数据，不调用系数函数。`previews/` 中的 PNG 是实际 Wolfram 图形导出；`top_angular_additional_2x2.png` 对应新增四图的组合。

## 验证

本次在 Wolfram 15.0.1 的新内核中实际计算并保存 notebook。原始计算约 145 秒；复用缓存时数据加载约 1 秒。

- Gauge：ν={0.01,0.1,0.6,1,2}，60 位比较，含数值精度估计的最大相对差异约 2.0×10^-60。
- Top：ν={0.1,0.2,0.73,1,2,5}，70 位比较，最大相对差异估计约 3.6×10^-67。
- 原第 8 节全部 40,401 个复数角网格样本：相对差异估计约 1.4×10^-30（原网格值保存至 30 位）。
- 实际幅度范围约 [4.6255592268,6.3827437009]×10^-6；主值相位约 [−2.3287056081,−2.1907038692] rad。
- 默认网格、有限复数值、相位分段、A3=−3A2、保存输入及新增 2×2 输出均已检查。
- 新增输入格已连同模块、参数格在另一个新内核执行通过，未加载原 notebook 的包，记录见 `notebook_replay.txt`。
- 十张图已逐张检查，坐标标题、图例、参数框及竖直色条完整。新 notebook 不展开完整数据表。
- 新增 π/2 图满足 x、y 独立反号对称，π 图与原 ϕ3=0 图沿 y 反射一致；包含数值精度估计的完整网格复数绝对残差低于 10^-28。
- 35 个原有源码/论文文件的 SHA256 全部保持不变，论文 PDF 未更新。主要参考文件哈希见 `source_hashes.json`。

详细数值及检查结果见 `validation_report.txt`，任意精度结果保存在 `validation_results.wxf`。这些检查针对移植与绘图复现。
