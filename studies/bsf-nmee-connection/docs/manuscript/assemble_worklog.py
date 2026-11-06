"""Create one current research log with archived KM and material records."""
from pathlib import Path
import hashlib
import json
import re

ROOT=Path(__file__).resolve().parent
REPO=Path('/Users/jiaqi/Desktop/worksofJiaqi/quantum-geometry-reproductions')
STUDY=REPO/'studies/bsf-nmee-connection'
MATERIAL=Path('/Users/jiaqi/Desktop/worksofJiaqi/project-x/outputs/MoSSe-NMEE')
archive=ROOT/'归档/研究记录'
entries=[]


def preserve(source, relative, label):
    source=Path(source)
    target=archive/relative
    target.parent.mkdir(parents=True,exist_ok=True)
    raw=source.read_bytes()
    text=raw.decode()
    # Keep access to underlying data after bringing a log into Obsidian.
    def link(match):
        name, href=match.groups()
        href=href.strip('<>')
        if '://' in href or href.startswith('#'):
            return match[0]
        base, _, anchor=href.partition('#')
        path=(source.parent/base).resolve()
        if path.exists():
            tail=('#'+anchor) if anchor else ''
            return f'[{name}](<{path}{tail}>)'
        return match[0]
    text=re.sub(r'\[([^\]]*)\]\(([^)]+)\)',link,text)
    header=f'> 历史研究记录，保留当时结论；当前判断以根目录 worklog 为准。[原始记录](<{source}>)。\n\n'
    target.write_text(header+text)
    entry={'label':label,'source':str(source),'archive':str(target.relative_to(ROOT)),
           'original_sha256':hashlib.sha256(raw).hexdigest(),
           'link_policy':'existing relative data links made absolute; original source unchanged'}
    entries.append(entry)
    return f'[{label}]({entry["archive"]})'


km=[]
for stage in range(2,8):
    km.append(preserve(STUDY/f'docs/STAGE{stage}_WORKLOG.md',f'KM/阶段{stage}.md',f'KM 阶段 {stage}'))
km8=preserve(STUDY/'docs/STAGE8_PHYSICS.md','KM/阶段8_极化与幂次.md','KM 极化与幂次')
km9=preserve(STUDY/'docs/STAGE9_MODEL_PHYSICS.md','KM/阶段9_完整路径.md','KM 完整滑移路径')
local=preserve(STUDY/'docs/FIXED_STACK_LAYER_CXX_WORKLOG.md','KM/固定堆垛的层响应.md','固定堆垛的层响应')
env=preserve(STUDY/'docs/DOMAIN_ENVIRONMENT_WORKLOG.md','KM/器件环境与投影零点.md','KM 环境响应')
initial=preserve(MATERIAL/'MoSSe-NMEE-结果.md','MoSSe/初始材料对照.md','MoSSe 初始对照')
distance=preserve(MATERIAL/'7269-mechanism/WORKLOG.md','MoSSe/S-S层距与带间分解.md','S/S 层距与带间分解')
orbital=preserve(MATERIAL/'7269-mechanism/orbital-model/WORKLOG.md','MoSSe/轨道模型与门控.md','轨道模型与门控')
blind=preserve(MATERIAL/'7266-blind/WORKLOG.md','MoSSe/S-Se堆垛的独立预测.md','S/Se 独立预测')
mirror=preserve(MATERIAL/'ss-mirror-domains/WORKLOG.md','MoSSe/镜像伙伴.md','MoSSe 镜像伙伴')
readout=preserve(MATERIAL/'domain-readout/WORKLOG.md','MoSSe/器件环境与实验读出.md','MoSSe 环境与读出')
story=preserve(STUDY/'docs/story-20261106/WORKLOG.md','稿件/最初物理叙事.md','最初物理叙事')

document=fr'''# worklog｜Layer NMEE

> 本项目唯一的当前工作记录。正文、模型定义和全部讨论图见 [正文](正文.md)。这里集中记录已经做过什么、哪些结论仍成立、哪些说法已更正，以及数据在哪里。历史记录的“下一步”只表示当时计划，不是当前待办。

## 当前物理状态

- 观察量是绝缘、固定离子条件下的静态 $E_xB_z$ 诱导层自旋；垂直控制电场或层势差与两个探测场分开定义。
- KM 与 MoSSe 均有两类现象：界面能够使指定层相对同朝向孤立层反向；总响应也能在两层自身不反向时过零。两者不能仅凭总信号的符号区分。
- 当前 KM 主例是 AB–BA 路径的 $\lambda=1/4$，邻近位置与小幅参数变化保留两类过程。KM AB/BA 作为高对称对照，不称它们已经重现同一门控补偿路径。
- MoSSe S/S 堆垛 A/B 由水平 $m_{{001}}$ 联系。它交换上下层并反转面内自旋；所以两个堆垛都可呈下负上正，而优势层相反。
- 尚未建立决定局域反向的可预测电子条件，也未证明实验可测性或 PRL 首创性。当前不扩材料，不自动启动新的远端作业。

## 2026-11-06｜连续稿件与目录整理

本轮按用户要求，只整理已有结果和文档，没有新增物理计算。

1. 新增正文图 2：将已有37个邻域计算点合并。上排比较 KM $\lambda=0.24,0.25,0.26$，下排比较 $m=44.9,45,45.1$ meV。每排依次展示下层随杂化翻向、门控时上下层保持异号、总响应过零。不同参数使用不同颜色和符号；总值面板放大过零区间，各层面板保留完整已算门控范围。
2. 原“正文”“补充讨论与图”“理论定义与公式核对”的物理内容并入同一篇正文，按问题—模型—邻域与对照—材料—实验含义—完整定义串联。所有图顺序编号1–10，保留图注。现在不划正文/SM，也不按PRL最终篇幅删减。
3. 日常入口固定为正文.md和worklog.md。图片放“附件”；原批注稿、分散说明、各轮稿件及旧图移入“归档/整理前-20261106”，保持原来的相对目录关系。没有删除研究内容。
4. 原来在计算目录里的KM与MoSSe记录已导入“归档/研究记录”，见下表。数据与计算输入保留在原计算目录；归档副本将能解析的旧数据链接指回原处。

## 必须继承的更正与约定

| 项目 | 当前应使用的解释 |
|---|---|
| KM AB/BA 与 MoSSe A/B | 前者是模型端点；后者是材料 S/S 镜像伙伴，不能混用名称或分数坐标。KM 四分之一构型的伙伴是四分之三构型。 |
| 旧 KM 约±10.25 meV 零点 | 原生 X 方向的纵向投影零点；横向分量仍非零。旧环境扫描的数值沿用该投影，不能改称完整自旋相消。 |
| 当前 KM 约4.7 meV 零点 | 电场沿保留竖直镜面的面内法线，磁场沿 z；另外两个自旋分量被禁止，因此是此 EB 阶的完整诱导自旋相消。旧点与新点属于不同探测方向。 |
| 两种镜面 | 竖直镜面保留给定构型、不交换层，限制自旋方向；水平 m001 联系两个伙伴、交换层。 |
| 层势差与栅场 | 模型 U 是内部层势能差，不能直接换算为材料栅压；材料 Fz 的正号沿用原计算电子受力场约定。 |
| 稳定性 | KM 电子 Hamiltonian 没有堆垛总能量。BiDB的滑移稳定不是所有声子均已验证稳定，更不是已完成写入/切换实验。 |
| 原NMEE公式 | 采用完整静态混合响应；原作者已确认的漏项不作为拟合目标。与驱动DC分支、磁矩归一化保持区分。 |
| 轨道归因 | 全部 d+p 层间块的干预比“唯一 S pz 键决定翻号”更稳健；单轨道删除结果依赖子空间定义。 |

## 研究沿革：从 KM 到材料和实验读出

| 阶段 | 已得到的物理、目前如何使用 | 完整记录 |
|---|---|---|
| KM 位移相关模型 | 从两张镜像相关单层出发建立层间跃迁；在绝缘条件下计算完整静态自旋响应。 | {km[0]} |
| KM 滑移与通道 | 得到有隙响应变号；不同层间跃迁通道不可简单相加。有限耦合的非加性不能全部解释为弱耦合干涉项。 | {km[1]}、{km[2]}、{km[3]} |
| KM 对称性与低能模型 | 参考构型附近的迹/方向差有不同最低阶；谷区模型可重现指定参数下的主要响应与反转。此为早期参数组，不与当前m=45的绝对数值混合。 | {km[4]}、{km[5]} |
| KM 电极化与完整路径 | 建立电子极化参照和完整滑移；同号电偶极不能单独决定自旋分量的方向。没有由电子模型证明铁电能垒或稳定堆垛。 | {km8}、{km9} |
| MoSSe 初始对照 | 孤立单层、反演双层、S/S双层给出不同的层响应与总值；反演参照的总值由对称性禁止。 | {initial} |
| S/S 层距控制 | 同朝向下层从接触时负响应恢复为正；较深与顶部占据带贡献竞争，不能只看Γ点VBM/CBM。 | {distance} |
| 材料轨道模型及门控 | 仅改变成组层间杂化可重现同层反向；原层距门控的已算区间则主要改变层间补偿。唯一最小轨道机制仍未确定。 | {orbital} |
| S/Se 新堆垛 | 封存下层负响应预测后，独立DFT支持其方向及明显增强。不同终端/取向不能简称单纯滑移。 | {blind} |
| 固定 KM 的两类反转 | 在m=45附近找到局域分量反转与总投影补偿；随后完整张量核查修正了“整个自旋”的表述。 | {local} |
| MoSSe 镜像伙伴 | 独立DFT给出两层交换并反号的关系，总响应零点约为±0.24 V/nm；上下层自身方向在该门控区间保持。 | {mirror} |
| 环境与仪器 | 公共标量偏置可平移曲线；真实单侧扰动还会改幅值与间距。探针权重可使仪器零点不同于总自旋零点；理想磁读出示例量级极小。 | {env}、{readout} |
| 稿件形成与批注 | 从分散结果转成模型—材料叙事，补Hamiltonian、带组标号、图注；此后改回四分之一主例并核清探测方向。 | {story}；[第一轮修订](归档/整理前-20261106/WORKLOG-20261106-批注修订.md)、[第二轮修订](归档/整理前-20261106/WORKLOG-20261106-第二轮修订.md)、[第三轮修订](归档/整理前-20261106/本轮工作记录-20261106c.md) |

## 同行阅读与用户批注

LEE、band-edge mixture、layer-dependent disorder、蓝磷界面、pentagonal-chain以及原NMEE正文/SM的阅读记录集中保留在[文献核查笔记](归档/整理前-20261106/文献图式核查.md)。这些记录区分文献做到了什么、我们借鉴什么、以及哪些机制尚未建立。它们不是本轮新增调研。

[原始带批注稿](归档/整理前-20261106/正文讨论稿.md)保持原文；[上次连续修订前的正文](归档/整理前-20261106/正文讨论稿-修订版.md)及配套说明也全部保留。当前阅读入口只使用本目录的正文.md，避免把旧版安排误作当前结论。

## 数据与接手位置

- [KM与手稿仓库](<{REPO}>)：基础模型、原始结果、历史核查；当前连续稿对应 studies/bsf-nmee-connection/docs/manuscript。
- [MoSSe计算记录](<{MATERIAL}>)：结构、模型干预、独立DFT汇总、镜像和读出结果。本轮没有移动原计算目录，也没有复制受许可的VASP程序、势文件或密钥。
- [本轮汇总图源数据](附件/km_neighbourhood.csv)：37行均来自此前已算的两个表；数值保持原样，source列标明来源。
- [本轮计算与整理核查](附件/整理核查.json)：图源一致性、图号/公式/链接、归档前后文件校验。此前精细四分之一计算与邻域检查使用不同积分精度，图2内部始终使用同一邻域精度。

## 下一轮

按用户要求，先完成稿件和目录整理，再在对话中讨论项目架构与物理新颖性；新的方向与算例随讨论确定。后续推进只更新这两个稳定入口；历史细节进入归档，计算资料留在计算仓库，不再在Obsidian根目录累加新的并列说明文件。
'''

(ROOT/'worklog.md').write_text(document)
(ROOT/'history_manifest.json').write_text(json.dumps(entries,indent=2,ensure_ascii=False)+'\n')
print('Unified worklog assembled; imported',len(entries),'historical records.')
