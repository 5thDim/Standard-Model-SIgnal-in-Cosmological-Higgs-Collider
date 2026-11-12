# Signal Plots

打开根目录的 `Signal Plots.nb` 即可查看十张已运行的图。正常运行只加载本目录中的五个模块，不读取原 notebook、原辅助脚本或论文。原文件均保留。

| 文件 | 职责 |
|---|---|
| `coefficients.wl` | 最终解析有限部、特殊函数稳定求值、复数系数接口及总角向组合 |
| `plot_data.wl` | 精确默认网格、复数扫描、固定质量扫角、幅度/相位投影、WXF 保存与缓存 |
| `plot_style.wl` | Latin Modern Roman 字体、字号、线宽、图尺寸、留白、图例、渐变色及分段色条 |
| `nu_scan.wl` | 根据已算好的复数扫描画 ln\|C\| 和主值相位；相位跳变处断线 |
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
