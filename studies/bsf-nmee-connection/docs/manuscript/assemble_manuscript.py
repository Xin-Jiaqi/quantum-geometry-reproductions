"""Consolidate the accepted physical text; keep review history in the worklog."""
from pathlib import Path
import re
import shutil
import json

ROOT = Path(__file__).resolve().parent
PREV = ROOT.parent/'revision-20261106c'
OUT = ROOT/'附件'
OUT.mkdir(exist_ok=True)
(ROOT/'data').mkdir(exist_ok=True)
for name in ['tensor_check.json','normal_main.csv','gate_extra_main.csv']:
    shutil.copy2(PREV/name,ROOT/'data'/name)
figure_names = {
    'Fig1_KM':'Fig01_KM_overview',
    'FigS4_KM_AB_BA':'Fig03_KM_AB_BA',
    'FigS6_probe_and_range':'Fig04_KM_probe_direction',
    'Fig2_MoSSe_distance':'Fig05_MoSSe_distance',
    'FigS5_orbital_levels':'Fig06_MoSSe_orbitals',
    'Fig3_MoSSe_controls':'Fig07_MoSSe_gate_controls',
    'FigS1_environment':'Fig08_MoSSe_environment',
    'FigS2_probe_weight':'Fig09_MoSSe_probe_weight',
    'FigS3_numerical_checks':'Fig10_MoSSe_numerical_checks',
}
for old,new in figure_names.items():
    for ext in ['png','pdf','svg']:
        shutil.copy2(PREV/'figures-20261106c'/f'{old}.{ext}',OUT/f'{new}.{ext}')

numbers={'1':'1','2':'5','3':'7','S1':'8','S2':'9','S3':'10','S4':'3','S5':'6','S6':'4'}
def remap(text):
    text=re.sub(r'图\s*(S\d+|\d+)',lambda m:'图 '+numbers.get(m[1],m[1]),text)
    text=re.sub(r'(?<![A-Za-z])S[1-6](?!\d)',lambda m:numbers[m[0]],text)
    for old,new in figure_names.items():
        text=text.replace(f'figures-20261106c/{old}.png',f'附件/{new}.png')
    return text

def section(text, title, next_title=None):
    begin=text.index(title)+len(title)
    end=text.index(next_title,begin) if next_title else len(text)
    return text[begin:end].strip()

main=(PREV/'正文讨论稿-修订版.md').read_text()
theory=(PREV/'理论定义与公式核对.md').read_text()
supp=(PREV/'补充讨论与图.md').read_text()
intro=section(main,'## 1．问题：界面改变一层的自旋方向，还是改变两层的相互抵消？','## 2．双层 Kane–Mele 模型')
model=remap(section(main,'## 2．双层 Kane–Mele 模型','## 3．MoSSe：检验界面诱导的层响应反转'))
material=remap(section(main,'## 3．MoSSe：检验界面诱导的层响应反转','## 4．物理意义与下一步机制问题'))
meaning=remap(section(main,'## 4．物理意义与下一步机制问题','## 参考文献与阅读路径'))
references=section(main,'## 参考文献与阅读路径').split('本稿的图文顺序')[0].strip()

model=model.replace('数值所用键方向及完全展开的 Bloch 矩阵见理论说明。','键方向及完全展开的 Bloch 矩阵列于第 6 节。')
model=model.replace('具体位移、坐标与层间跃迁形式在理论说明中列全。','具体位移、坐标与层间跃迁形式列于第 6 节。')
model=model.replace('旧稿使用的另一电场方向及其非零横向响应保留在图 4，避免混用两套系数。','一般电场方向下，纵向投影过零仍可能留下横向自旋；图 4 对比这两种探测条件。')
model=model.replace('与原 NMEE 几何表达式的对应及对角项补全放入理论说明；','与原 NMEE 几何表达式的对应及完整静态导数见第 6 节；')
model=model.replace('［图 1(d)及补充表］','［图 2(a–c)］')
model=model.replace('小幅改变 $m$ 也不立即消除它们。','将 $m$ 小幅改为 44.9、45.1 meV 也保留两类过程［图 2(d–f)］。')
model=re.sub(r'四分之一位置没有额外的对称性特殊性；.*?这里的门控效应同样需要层间作用。',
             lambda _:r'为检验上述对照是否依赖精确的四分之一位置，我们在同一探测方向下比较 $\lambda=0.24,0.25,0.26$，并在固定 $\lambda=0.25$ 时比较 $m=44.9,45,45.1$ meV。',model)

neighbour_figure=r'''![图 2](附件/Fig02_KM_neighbourhood.png)

图 2｜两种反转过程在邻近 KM 堆垛与参数中共存。上排固定 $m=45$ meV，颜色与符号区分 $\lambda=0.24,0.25,0.26$；下排固定 $\lambda=0.25$，区分 $m=44.9,45,45.1$ meV。(a,d) 固定 $U=0$，增强层间杂化使同一下层由负转正。(b,e) 固定 $t_\perp=10$ meV，改变层势差：实线为下层，虚线为上层，各层保持原方向。(c,f) 同一组数据的总响应，放大 $U=4$–$5$ meV 的过零区间。系数单位统一为 $10^{-5}(\mathrm{nm\,V\,T})^{-1}$；连线仅连接计算点，各参数的取样数量不同，未拟合补点。

图 2 将“层自身翻向”和“层间补偿”分别显示出来：三种相邻堆垛都保留前者，也都能在上下层保持异号时实现后者。变化的是临界控制量与幅值，两个过程并不依赖精确的四分之一位置。这里验证的是邻域中的共存，而非整个滑移路径的普遍定律；更远的 $\lambda=0.30$ 或 $m=46$ meV 未在相同杂化扫描区间重现下层反向。这些对照也说明，对称性允许一种响应，并不保证任意电子参数下都能把它调过零。
'''
anchor='水平镜像还给出可迁移的结构关系：'
model=model.replace(anchor,neighbour_figure+'\n'+anchor)
ab=remap(section(supp,'## S4．KM AB/BA：局域反向不要求四分之一构型','## S5．轨道投影给出了哪些信息？'))
probe=remap(section(supp,'## S6．投影零点与完整自旋相消必须区分','## S7．四分之一选择不是单点巧合，也尚非普遍定律'))
probe=probe.replace('旧的原生 $X$ 电场方向','沿模型原生 $X$ 的电场方向').replace('新主图的镜面法线探测方向','沿竖直镜面法线的电场方向')
probe=probe.replace('新的主图','图 1').replace('新主图','图 1')
probe=probe.replace('两种探测方向的数据分别计算，不互换系数名称。','两种电场方向分别计算。')
probe=probe.replace('旧的约 10.25 meV 过零','沿原生 $X$ 的约 10.25 meV 过零')
probe=probe.replace('详细变换见理论说明第 9 节。','坐标与轴矢量变换见第 6 节。')
model+='\n\n### 高对称构型与探测方向的对照\n\n'+ab+'\n\n'+probe
model=model.replace('图 1利用','图 1 利用')
model=model.replace('更大 $U$ 下，单层自身也会反向。','更大 $U$ 下，单层自身也会反向。对完全解耦且占据不变的绝缘层，刚性层势差不改变响应；这里的门控作用也依赖层间耦合。')

orbital=remap(section(supp,'## S5．轨道投影给出了哪些信息？','## 带组的定义和解释边界'))
orbital=orbital.replace('不是新的零磁场连续能带计算。','其物理对象是准备场下的离散轨道投影，而非零场连续能带。')
orbital=orbital.replace('目前可画的是离散能级／轨道投影，而非连续路径 fatband。连续路径、零场轨道投影和主导响应子空间应在下一轮机制计算中对应起来；不能仅为得到类似文献的图形而在离散点间补造色散。','轨道投影确定了候选态的成分，但未给出决定响应符号的矩阵元。连续能带上的轨道成分还需与主导响应子空间对应，才能据此构建更小的有效模型。')
orbital=orbital.replace('\n\n轨道投影确定了候选态的成分，但未给出决定响应符号的矩阵元。连续能带上的轨道成分还需与主导响应子空间对应，才能据此构建更小的有效模型。','')
material=material.replace('![图 7]',orbital+'\n\n![图 7]',1)
material=material.replace('非对称环境及实验读出的限制单列在[补充讨论与图](补充讨论与图.md)，不在当前结果段落中将这一关系直接宣称为已实现的畴读出方案。','这些过零场能否用于实际畴读出，还需检查环境与探针的空间权重，见第 4 节。')

environment=remap(section(supp,'## S1．单侧环境扰动影响过零场的绝对位置','## S2．磁探针的空间权重会改变信号的零点'))
weight=remap(section(supp,'## S2．磁探针的空间权重会改变信号的零点','## S3．过零场与端点符号的数值检查'))
numerics=remap(section(supp,'## S3．过零场与端点符号的数值检查','## 仍有价值、但暂不合并为主图的结果'))
groups=remap(section(supp,'## 带组的定义和解释边界','## S6．投影零点与完整自旋相消必须区分'))

def th(start, end=None):
    s=section(theory,start,end)
    return remap(s).replace('本轮','此处').replace('修订正文','正文').replace('本稿主图','图 1')

bloch=th('## 2．单层 Bloch 矩阵的完全展开','## 3．层间矩阵、平移及电子数')
bloch=bloch.replace('见第 9 节','见本节后文')
tunneling=th('## 3．层间矩阵、平移及电子数','## 4．静态系数的推导与机制分解')
tunneling=tunneling.replace('正文用“固定堆垛构型，保持半填充与绝缘态”，替换“固定配准和填充”。','')
response=th('## 4．静态系数的推导与机制分解','## 5．本轮复核范围')
paper=th('## 1．原 NMEE 的定义与当前系数','## 2．单层 Bloch 矩阵的完全展开')
paper=paper.replace('第一作者已在用户提供的回复中确认：','原 NMEE 第一作者曾确认：')
paper=paper.replace('归入推导背景，不作为本文机制主线；','属于响应定义的背景；')
horizontal=th('## 6．KM AB/BA 及中间构型的水平镜像','## 7．不依赖整体符号约定的参照')
horizontal=horizontal.replace('此处没有新的 DFT。','')
vertical=section(theory,'KM AB–BA 连线各构型保留竖直镜面。','## 10．四分之一构型的选择与检验范围')
vertical='KM AB–BA 连线各构型保留竖直镜面。'+vertical
vertical=vertical.replace('本稿主图','图 1').replace('不是把旧 $C_{XX}$ 曲线换一个标签。','不同方向的系数必须同时变换读出自旋和电场顶角。')
vertical=vertical.replace('本轮也在数值中观察到，故未将那个方向用作完整补偿的主图。','该允许分量也出现在计算中，因而不能沿那个方向只以纵向过零判定完整补偿。')
vertical=vertical.replace('新坐标','上述坐标').replace('新过零邻域','过零邻域')
decoupled=th('## 8．解耦参照下的势差控制','## 9．为什么重新选择探测方向？')
formula_checks=th('## 5．本轮复核范围','## 6．KM AB/BA 及中间构型的水平镜像')
formula_checks=formula_checks.replace('复核 MoSSe 分解代码后，应修正上一版对空态截断的表述：','MoSSe 的两种分解定义不同：')
formula_checks=formula_checks.replace('这只核对公式与实现，不是新的材料收敛证明。','这些检查核对公式与实现；材料精度由下面的独立对照判断。')
formula_checks=formula_checks.replace('原文漏项及推导对应保留在本说明，不在正文反复讨论术语争议。','')

document='''# 双层绝缘体中层分辨自旋响应的界面调控

> 连续讨论稿，更新于 2026-11-06。当前将物理问题、模型、材料、理论定义和验证放在同一篇文档中；待科学故事确定后再压缩篇幅、决定哪些内容移入 SM。模型使用“KM AB/BA”，材料使用“MoSSe S/S 堆垛 A/B”。研究沿革与数据来源统一见 [worklog](worklog.md)。

## 1．物理问题：一层自身翻向，还是两层补偿改变？

'''+intro+'\n\n## 2．KM 模型：同一固定堆垛中的两种反转\n\n'+model+'''

## 3．MoSSe：从模型命题到材料验证

'''+material+'''

## 4．从响应过零到实验辨认：环境与探针的作用

模型与材料中的总自旋过零，是一个明确的响应条件。把它用于实验辨认堆垛，还需要知道外部环境怎样改变响应，以及仪器是否给上下层相同的权重。

### 内建偏场与单侧环境

'''+environment+'\n\n### 测量信号是否与总自旋具有同一个零点？\n\n'+weight+'''

## 5．已确立的物理与尚待解决的机制

'''+meaning+'''

## 6．模型、响应与验证的完整定义

以下保留判断上述物理所需的公式、算符和数值边界。它们与结果属于同一份稿件，目前不单独划作补充材料。

### 单层矩阵与层间跃迁

'''+bloch+'\n\n'+tunneling+'''

### 静态混合响应及其与 NMEE 量子几何的关系

'''+response+'\n\n'+decoupled+'\n\n'+paper+'''

### 两种镜面分别约束什么？

'''+horizontal+'\n\n'+vertical+'''

### 公式核对、材料精度与带组分解

'''+formula_checks+'\n\n'+numerics+'\n\n'+groups+'''

## 参考文献

'''+references+'\n'

# Standalone old-document references must not survive the consolidation.
document=document.replace('理论说明','第 6 节').replace('补充表','图 2')
document=document.replace('所考察的 AB–BA','所考察的 KM AB–BA')
document=document.replace('图 1采用','图 1 采用').replace('AB 和 BA 的','KM AB 和 KM BA 的').replace('AB 的较细计算','KM AB 的较细计算')
document=document.replace('未在相同杂化扫描区间重现下层反向','未在已检查的杂化范围重现下层反向')
document=document.replace('### 公式核对、材料精度与带组分解\n\n',
    '### 公式核对、材料精度与带组分解\n\n图 1 的四分之一参考曲线采用 81×80×120 的谷区积分参数和 0.03125 T 导数步长；图 2 的所有邻域对比统一采用 41×45×72 与 0.0625 T。两图中的中心参数因此不是同一次积分，不据幅值末位比较差异。三项积分参数不代表三维均匀 k 网格。\n\n')
document=re.sub(r'\n{4,}','\n\n\n',document)
(ROOT/'正文.md').write_text(document)
(ROOT/'manuscript_map.json').write_text(json.dumps({
    'source_revision':'revision-20261106c', 'figure_renaming':figure_names,
    'new_figure':'Fig02_KM_neighbourhood',
    'sections':['physical question','KM and robustness','MoSSe','experimental interpretation','mechanism status','complete definitions'],
    'omissions':'duplicated definitions and reviewer conversation removed from prose; source text preserved in archive/worklog',
    'new_calculation':False,
},indent=2,ensure_ascii=False)+'\n')
print('Continuous manuscript assembled; historical sources unchanged.')
