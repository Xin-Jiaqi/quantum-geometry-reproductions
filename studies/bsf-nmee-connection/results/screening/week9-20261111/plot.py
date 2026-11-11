"""Week 9: signed zero-gate spin directions and absolute changes on one scale."""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.lines import Line2D

here = Path(__file__).resolve().parent
source = here.parent / "review7268-20261111/matched_kernel_interfaces.json"
rows = [r for r in json.loads(source.read_text()) if r["bands"] == 112 and r["mesh"] == 6]
assert [r["interface"] for r in rows] == ["S/S", "S/Se"]
order = [1, 2, 0]  # lower, upper, entire bilayer
zero = np.array([r["Czero"] for r in rows])[:, order]
minus = np.array([r["Cminus"] for r in rows])[:, order]
plus = np.array([r["Cplus"] for r in rows])[:, order]
delta = abs(plus) - abs(minus)
assert np.allclose(plus[:, :2].sum(axis=1), plus[:, 2], atol=1e-7)
assert (np.sign(minus) == np.sign(zero)).all() and (np.sign(plus) == np.sign(zero)).all()
font_path = Path("/System/Library/Fonts/STHeiti Light.ttc")
if font_path.exists():
    font_manager.fontManager.addfont(str(font_path))
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Heiti TC", "Noto Sans CJK SC", "DejaVu Sans"],
                     "svg.fonttype": "none", "pdf.fonttype": 42,
                     "font.size": 10, "axes.unicode_minus": False, "axes.linewidth": .7})
colors = ["#20678b", "#bd6a28", "#343d49"]
names = ["下层", "上层", "整片"]
fig = plt.figure(figsize=(9.8,4.4), facecolor="white")
fig.text(.055,.94,"垂直电场使两层自旋都减弱，整片响应却可以增强",fontsize=17,weight="bold")
fig.text(.055,.865,"Janus MoSSe 双层   ·   面内电场沿 +x，磁场沿 +z，+z 由下层指向上层",color="#435267",fontsize=10)
a = fig.add_axes([.055,.245,.355,.46])
b = fig.add_axes([.52,.245,.44,.46])
a.set_title("a   零垂直电场下的自旋方向",loc="left",pad=19,weight="bold",fontsize=11)
b.set_title("b   增加垂直电场后的响应幅度变化",loc="left",pad=19,weight="bold",fontsize=11)
a.set(xlim=(-1.1,2.45),ylim=(-.55,1.82)); a.axis("off")
for k, label in enumerate(names):
    a.text(k,1.52,label,ha="center",color=colors[k],weight="bold")
for i, row in enumerate(rows):
    y = .93 - i
    a.text(-1.07,y,row["interface"],va="center",weight="bold")
    for k, val in enumerate(zero[i]):
        sign=np.sign(val)
        a.annotate("", xy=(k+.22*sign,y+.07), xytext=(k-.22*sign,y+.07),
                   arrowprops={"arrowstyle":"-|>","color":colors[k],"lw":1.8,"mutation_scale":14})
        a.text(k,y-.18,"+x" if val>0 else "−x",ha="center",va="center",color=colors[k],fontsize=9)
for i, row in enumerate(rows):
    for k, value in enumerate(delta[i]):
        y=(1-i)*1.05+(.23-k*.23)
        b.barh(y,value,height=.15,color=colors[k])
        # Keep very small changes visible using their numeric labels.
        b.text(value + (.028 if value>=0 else -.028),y,f"{value:+.2f}".replace("-","−"),
               ha="left" if value>=0 else "right",va="center",fontsize=8,color=colors[k])
b.axvline(0,color="#768299",lw=.7)
b.set(xlim=(-1.62,.52),ylim=(-.4,1.48),yticks=[1.05,0],yticklabels=["S/S","S/Se"],
      xticks=[-1.5,-1.,-.5,0,.5],xlabel="幅度变化  [$10^{-9}$ (nm·V·T)$^{-1}$]")
b.tick_params(axis="y",length=0,pad=5)
b.spines[["left","right","top"]].set_visible(False)
b.legend(handles=[Line2D([],[],color=c,lw=5,label=n) for c,n in zip(colors,names)],
         loc="upper right",bbox_to_anchor=(1.01,1.16),frameon=False,ncol=3,fontsize=8,
         handlelength=.8,columnspacing=.9)
fig.text(.055,.092,"左图箭头只表示方向。右图比较垂直电场从 −0.10 到 +0.10 V/nm 时，各层及整片的响应幅度变化。",fontsize=9,color="#435267")
fig.text(.055,.042,"固定离子，各场点独立自洽；两种界面均取 6 × 6 网格、112 带。",fontsize=9,color="#435267")
for name in ["MoSSe-week9.png", "MoSSe-week9.pdf", "MoSSe-week9.svg"]:
    fig.savefig(here/name, dpi=600, bbox_inches="tight", pad_inches=.15)
