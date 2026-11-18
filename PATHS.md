# 有限修正路径与组分限值

本文件描述待新增能力；V1 正常业务、数字规则、修订流程及失败钩子见 README.md。
现有 direction 只描述零点右导数，不是有限修正路径。
本题不扩展 archive.export；原方向导出和全部原接口继续工作。

## 调用入口

```python
from spectral.response.requests import trace
from spectral.response.limits import acceptable
from spectral.archive.reading import read

path_id = trace(c, parent_id, [-1, 1], 0, 1)
report_id = acceptable(c, path_id, {"k1": (None, "3/4")})
path = read(c, path_id)
report = read(c, report_id)
```

- trace(c, parent_id, d, lo, hi, fail_at=None) 返回 path 成果 ID。
  parent_id 必须为已发布 unmix 成果；同一 d 作用于其全部记录。
  d 为沿父通道顺序的 list/tuple；长度匹配通道，数值沿用 parse_rational。
  lo、hi 也沿用该精确数值规则，要求 0 <= lo <= hi；区间为闭区间。
- acceptable(c, path_id, bounds, fail_at=None) 返回 acceptable 成果 ID。
  bounds 是 dict，键是父字典组分 ID，值是长度为 2 的 list/tuple：
  (lower, upper)。非 None 端点沿用精确数值规则；None 表示该侧无限制。
  下界不能大于上界；未知组分、非法容器或端点整次拒绝。
  限值针对声明含量单位的数值，不换算单位。
- 未知成果 ID 抛 KeyError；已存在但种类错误的成果、非法业务输入抛 ValueError。
  数值拒绝 float、bool、None（只有限值端点的 None 例外）及非有限/非法字符串。
  文件和数据库错误保持原类型。
- fail_at 沿用 None、after_payload、after_index；未知阶段抛 ValueError。
  新成果负载和目录同事务发布；失败不留下新成果、孤立负载或成功回执。
  同一父成果和解析后的同一请求稳定重试返回同一 ID；限值键的输入顺序不影响身份。
  不增加进程崩溃恢复或多工作者要求。

trace 的路径计算与发布归响应业务模块；limits 模块负责限值请求及已保存路径的分析与发布，二者通过 path_id 交接。

## 完整路径

对每条父测量，在整个 [lo,hi] 上求 y(t)=y+t*d 的非负系数。
首先全局最小化 ||A*x-y(t)||^2，再在全部残差最优点中最小化 ||x||^2。
矩阵、坐标、单位及记录顺序固定来自指定父成果，不追随源当前指针。
沿用 V1 的全部维度和字典支持，包括零列、重复列、比例列及秩亏。
允许全量面枚举、全量重算及逐记录处理，没有新增性能目标。

read(c,path_id) 的 kind 为 "path"，parents 恰为
[{"type":"artifact","id":parent_id}]。data 至少有：
- parent_id、cohort_id；
- channels、channel_units、components、component_units，沿用父顺序；
- direction（解析后的 d）、lo、hi；
- rows，按父记录顺序；每行含 sample_id、measurement_revision、segments。

每个 segment 恰含 lo、hi、constant、slope；向量沿 components 顺序，
在该段内 x(t)=constant+t*slope，使用全局参数 t，不是 t-segment.lo。
上述所有数字读取时均为 Fraction，向量与 rows/segments 为 list。
当 lo < hi 时，各段长度为正，端点升序，无缝覆盖整个父区间；
相邻段在共享端点给出同一点值，且相同 constant/slope 的相邻段必须合并。
段表示的是仿射含量变化，不以内部面标签的变化强行分段。
当 lo == hi 时，每行只有一个 lo==hi 的段：
constant 为该点规范含量，slope 为全零向量。

保存的是完整连续分段仿射路径，而非待求解请求。
不能用零点切线、有限网格或保守 right_limit 当成全部路径。
无需新增重构/残差路径或路径导出格式。

## 限值消费

acceptable 必须从已保存的指定路径取得全部系数段、组分坐标及参数范围；
不得重新解混、重算方向或用当前源替代它。发布报告后读取同样不计算。
每行同时满足 bounds 中全部上下闭限值的参数集合，按升序保存为全部最大闭区间。
相接区间合并，孤立合格点表示为 [t,t]；不取外包区间，不遗漏单点。
无合格参数为 []；空 bounds 返回 [[lo,hi]]，包括 lo==hi 情形。

read(c,report_id) 的 kind 为 "acceptable"，parents 恰为
[{"type":"artifact","id":path_id}]。data 至少有：
- path_id；
- bounds，保存该请求的组分键及 [lower,upper]，有限端点解码为 Fraction；
- rows，顺序同路径，每行含 sample_id、measurement_revision、intervals。
intervals 为由两个 Fraction 组成的 list 的 list。

路径及限值成果均保存独立输入和明确父链。
字典/测量更正、调用方输入修改、读取对象修改、关闭重开不改变旧成果；
read 只解码，限值分析不因重新求解而掩盖不完整路径。
原 SQLite 表、修订、逻辑身份和 V1 旧负载保持兼容。

## 手算及运行

path_examples.py 提供固定输入、调用及预期数据，没有路径或限值算法。
finite_boundary 的首条测量 A=[[1,1],[0,1]]、y=[1,0]、d=[-1,1]：
[0,1/2] 上 x=[1-2t,t]；[1/2,1] 上 x=[0,1/2]。
k1<=3/4 接受 [0,1]；k1==1/4 仅接受 [1/4,1/4]；
k1==1/2 接受 [1/2,1]。第二条测量 y=[0,1] 全段为 [0,1/2]。
duplicate_simultaneous 核对重复列同时归零；其余样例核对比例列、零矩阵、
非零左端及单点父区间。样例不是完整验收，也不限制其他合法输入。

基线命令（从项目根运行）：
```
python3 -B verify.py
python3 -B demo.py
python3 -B path_examples.py
```
最后一条只核对手算输入与已有单点能力，明确不表示 V2 已实现。
完成新增能力后可运行 python3 -B path_examples.py --v2，查看真实路径及限值结果。

交付时沿现有 .pointwise-runtime.json 联合接通全部原回归与新增独立核账，
demo 保留原业务演示并展示跨断点、孤立点、更正后重开旧成果。
清单保留 version=1、kind=cli；setup 为 argv 数组列表，verify/demo 各为 argv 数组。
独立核账应检查分段及断点的分层最优性和完整覆盖，并有简短的完整性依据，
不能仅用少量参数采样代替整段核验。还须核对多记录身份顺序、非法输入、
两种发布失败阶段、稳定重试、输入/返回值隔离、旧成果重开及读取不求解。
复用原失败钩子与回归设施，保留原测试字节；交接说明解释路径生产与限值消费
如何共同保持坐标、全局参数及闭端点语义。
