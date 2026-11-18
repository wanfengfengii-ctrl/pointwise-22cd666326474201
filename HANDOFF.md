# 交接说明：有限修正路径与组分限值

本文件说明新增的 `trace` / `acceptable` 能力如何与 V1 共存，以及路径生产与限值
消费如何共同保持坐标、全局参数和闭端点语义。运行方式见 README.md、PATHS.md 与
`.pointwise-runtime.json`。

## 模块与调用链

| 模块 | 职责 |
| --- | --- |
| `spectral/response/paths.py` | 从父 unmix 负载生产完整分段仿射路径并发布 `path` 成果 |
| `spectral/response/limits.py` | 解析限值请求、只消费已保存路径、发布 `acceptable` 成果 |
| `spectral/response/requests.py` | 公开入口 `trace`、`acceptable`，负责数值与容器校验 |

```python
from spectral.response.requests import trace
from spectral.response.limits import acceptable
from spectral.archive.reading import read

path_id = trace(c, parent_id, [-1, 1], 0, 1)
report_id = acceptable(c, path_id, {"k1": (None, "3/4")})
```

`trace` 在 `requests.py` 解析输入、校验维度与区间，调用 `paths.analyze_path`
逐记录生产路径，再由 `paths.publish_path_result` 发布。`acceptable` 在
`limits.py` 内部由 `read` 取回 `path` 成果，交给 `parse_bounds` 校验限值，
再由 `build_acceptable_payload` 逐行计算并发布。两段通过 `path_id` 交接，
限值分析不再触碰矩阵、测量或求解器。

## 路径生产：为什么是完整路径

对每条父测量，在整个闭区间 `[lo,hi]` 上令 `y(t) = y + t*d`，先全局最小化
`||A x − y(t)||²`，再在全部非负残差最优点中最小化 `||x||²`。

- 对每个支撑集 `S`，把最优解一次解成 `x(t) = constant + t*slope` 的仿射族
  （对 Gram 矩阵做最小范数求解；同一组方程同时解 `y` 与 `d` 的右端项）。
- 该仿射族的非负约束 `constant_j + t*slope_j ≥ 0` 给出可行的有理参数区间
  `[feasible_lo, feasible_hi]`。
- 规范解对 `t` 唯一且连续（目标函数严格凸、可行集非空闭凸）。在固定支撑上它是
  有理系数的仿射函数，因此公式只可能在两类有理点改变：
  1. 某支撑某分量的系数过零，`t = −constant_j / slope_j`；
  2. 两个支撑的仿射表达式重合，`constant_a + t*slope_a ≡ constant_b + t*slope_b`
     的解（唯一候选 `t`）。
- 这些点全部并入断点集合，再加上 `lo`、`hi`，切成有限个片；片内不会出现隐藏的
  公式切换。
- 每个片用**整段解析证明**选胜者并复核：对每个与片有开区间交叠的支撑，在其自身
  可行交叠上证明残差平方和之差在整个闭片非负（有理二次函数，端点加顶点即精确
  符号），残差恒等时再证明范数平方差非负；残差差的有理零点处补验范数，处理并列。
  不是采样、不是零点切线、不是保守 `right_limit`。
- 断点只由公式变化决定，不因内部面标签变化强行分段；相邻相同 `(constant, slope)`
  合并。`lo == hi` 时直接调用 `solve_canonical` 得单点规范含量，`slope` 置零。

`canonical_path` 对 `lo > hi` 抛 `ValueError`；`lo < hi` 时保证每段长度为正、
端点升序、无缝覆盖、相邻共享端点同值、相同式必合并。

## 限值消费：坐标与闭端点

`acceptable` 从路径成果读取 `components` 顺序与每行 `segments`。对每个限值键用
`components` 中的位置索引 `constant` 与 `slope`（绝不按支撑顺序或键序猜坐标），
把每个段裁剪到同时满足全部上下闭限值的参数子区间：

- 每个仿射约束 `lower ≤ constant_j + t*slope_j ≤ upper` 精确变换为一个闭
  `t` 区间；`slope_j == 0` 时退化为常量判定。
- 无界侧不裁剪；空 `bounds` 直接保留每个段原区间，`lo == hi` 时得 `[[lo,lo]]`。
- 段裁剪结果升序排序并合并相接/重叠区间；无合格参数为 `[]`，孤立合格点保留为
  `[t,t]`，不取外包区间、不漏单点。
- `bounds` 发布为按键排序的 `{组分: [lower, upper]}`，有限端点为 `Fraction`，
  `None` 保持 `None`；逻辑身份基于排序后的键值，故键的输入顺序不影响 `ID`。

所有端点均为 `Fraction`，`read` 返回的 `intervals` 为 `[[Fraction, Fraction], …]`。

## 成果、身份与失败语义

- `path` 成果 `parents` 恰为 `[{"type":"artifact","id":parent_id}]`，`data` 含
  `parent_id`、`cohort_id`、`channels`、`channel_units`、`components`、
  `component_units`、`direction`、`lo`、`hi`、`rows`（按父顺序，含
  `sample_id`、`measurement_revision`、`segments`）。
- `acceptable` 成果 `parents` 恰为 `[{"type":"artifact","id":path_id}]`，`data`
  含 `path_id`、`bounds`、`rows`（含 `sample_id`、`measurement_revision`、
  `intervals`）。
- 路径只来自父成果保存的 `input` 快照，不读取字典/测量当前指针；`read` 仅解码，
  限值分析只用已存段。更正、对象修改、关闭重开都不改变旧成果；同一父和同一解析
  请求稳定重试返回同一 `ID`。
- 失败沿用原设施：`publish_artifact` 的 `after_payload` / `after_index` 抛
  `RuntimeError` 并回滚，未知阶段 `ValueError`。未知成果 `KeyError`，种类错误或
  非法业务值 `ValueError`。数值沿用 `parse_rational`：拒绝 `float`、`bool`、
  `None`（仅限值端点的 `None` 例外）、非有限/非法字符串。
- 未扩展 `archive.export`；原方向导出、原接口和原测试字节全部保留。

## 坐标、全局参数、闭端点的一致性依据

1. **坐标**：路径向量与 `components` 同序，限值用 `components` 名字定位索引；
   两者都取自同一父成果，没有隐式重排。
2. **全局参数**：每段 `x(t) = constant + t*slope` 用全局 `t`，不是
   `t − segment.lo`；限值裁剪与断点合并都在同一全局坐标上。
3. **闭端点**：断点显式加入区间两端，段端点升序且共享端点同值；限值裁剪保留
   `≤`/`≥` 的闭边界，因此孤立合格点能作为 `[t,t]` 出现而不被丢弃。

## 核账

`tests/test_paths.py` 独立于生产代码：用行列式/Cramer 自行解各支撑的仿射族，
在整段上用有理二次精确符号证明分层最优性、完整覆盖、断点与合并；并覆盖多记录
身份顺序、非法输入、两种发布失败阶段、稳定重试、输入/返回值隔离、旧成果重开、
读取不求解（用 `mock` 断言 `canonical_path` 与 `solve_canonical` 未被调用）。

`path_examples.py` 保留固定手算样例；`--v2` 用真实 API 复核 `finite_boundary`、
`duplicate_simultaneous`、`proportional`、`zero_matrix`、`singleton`、
`nonzero_left`。`demo.py` 在原业务演示外展示跨断点、孤立点、更正后重开旧成果。
