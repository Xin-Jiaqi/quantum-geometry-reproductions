"""Compare zero-gate layer NMEE of matched isolated MoSSe sheets and their S/Se bilayer."""
from pathlib import Path
import hashlib
import json

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import MultipleLocator
import numpy as np

here = Path(__file__).resolve().parent
states = here.parent / "review7268-20261111/states"
structures = here.parent / "overnight7268-20261110/compute_inputs"


def response(mesh, system, bands):
    state = json.loads((states / f"k{mesh}_{system}0/result.json").read_text())
    structure_hash = hashlib.sha256((structures / f"POSCAR.{system}").read_bytes()).hexdigest()
    assert state["Fz_V_nm"] == 0 and state["complete"] and state["input_sha256"]["POSCAR"] == structure_hash
    return np.array(state["responses"][str(bands)]["standard"]["C"])[:, 0, 0]


def comparison(mesh):
    lower = response(mesh, "lower", 64)[0]
    upper = response(mesh, "upper", 64)[0]
    isolated = np.array([lower, upper, lower + upper])
    coupled = response(mesh, "A", 112)[[1, 2, 0]]
    np.testing.assert_allclose(coupled[2], coupled[:2].sum(), atol=1e-7)
    assert isolated[0] > 0 > coupled[0] and isolated[1] > 0 and coupled[1] > 0
    assert np.all(abs(coupled[:2]) > abs(isolated[:2])) and 0 < coupled[2] < isolated[2]
    return isolated, coupled


for mesh in (6, 9):
    comparison(mesh)
isolated, coupled = comparison(12)
font_manager.findfont(font_manager.FontProperties(family="Times New Roman"), fallback_to_default=False)
plt.rcParams.update({"font.family": "Times New Roman", "font.size": 10,
                     "pdf.fonttype": 42, "svg.fonttype": "none",
                     "axes.spines.top": False, "axes.spines.right": False,
                     "mathtext.fontset": "custom", "mathtext.rm": "Times New Roman",
                     "mathtext.it": "Times New Roman:italic", "mathtext.bf": "Times New Roman:bold",
                     "mathtext.cal": "Times New Roman:italic", "mathtext.fallback": None,
                     "axes.linewidth": .8, "axes.unicode_minus": True})
fig, ax = plt.subplots(figsize=(3.95, 3.25))
fig.subplots_adjust(left=.18, bottom=.15, right=.98, top=.85)
positions = np.arange(3)
width = .28
for offset, values, color, edge, label in (
        (-.17, isolated, "#dededb", "#757773", "Isolated sheets"),
        (.17, coupled, "#2b6485", "#225774", "S/Se bilayer")):
    bars = ax.bar(positions + offset, values, width, color=color, edgecolor=edge,
                  linewidth=.8, label=label, zorder=2)
    for bar, value in zip(bars, values):
        ax.annotate(f"{value:+.1f}".replace("-", "−"),
                    (bar.get_x()+bar.get_width()/2, value),
                    xytext=(0, 3 if value > 0 else -3), textcoords="offset points",
                    ha="center", va="bottom" if value > 0 else "top", fontsize=9)
ax.axhline(0, color="#484b4c", linewidth=.75, zorder=1)
ax.set(xticks=positions, xticklabels=["Lower layer", "Upper layer", "Net"],
       xlim=(-.55, 2.55), ylim=(-72, 83),
       ylabel="Spin response, $C_{xx}$  [$10^{-9}$ (nm·V·T)$^{-1}$]")
ax.yaxis.set_major_locator(MultipleLocator(20))
ax.tick_params(direction="out", length=3)
ax.legend(loc="lower center", bbox_to_anchor=(.5,1.01), ncol=2, frameon=False,
          handlelength=1.15, handletextpad=.5, columnspacing=1.3)
for name in ("MoSSe-layer-response.pdf", "MoSSe-layer-response.png", "MoSSe-layer-response.svg"):
    path = here / name
    fig.savefig(path, dpi=600)
    if path.suffix == ".svg":
        path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines())+"\n")
