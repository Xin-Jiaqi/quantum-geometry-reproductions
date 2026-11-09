#!/usr/bin/env python3
"""Plot the frozen computed angle scan; recompute with reproduce_core.py first."""
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
r=json.loads((ROOT/'results/core.json').read_text())['angle_scan']['data']
x=np.degrees([a['phi_rad'] for a in r])
plt.rcParams['svg.fonttype']='none'
fig,ax=plt.subplots(figsize=(7,3.6))
for protocol,label in [('printed_D10','Orbital: printed D10'),('density_A9','Orbital: density A9')]:
    ax.plot(x,[a[protocol]['orbital'][0]*1000 for a in r],marker='o',markersize=3,label=label)
ax.plot(x,[a['density_A9']['spin'][0]*1000 for a in r],linestyle='--',label='Spin: both protocols')
ax.set(xlabel='Néel angle (degrees)',ylabel='Im χ(x,z) [mA / (m T)]',
       title='CuMnAs model: magnetic-order reversal',xticks=[0,90,180,270,360])
ax.legend(fontsize=8);ax.grid(alpha=.25)
fig.text(.5,.01,'μ = 0.08 eV; T = 25 K; ℏω = 0.1 meV; a = 3.82 Å (explicit input); +ω harmonic',ha='center',fontsize=8)
fig.tight_layout(rect=[0,.05,1,1])
(ROOT/'figures').mkdir(exist_ok=True)
fig.savefig(ROOT/'figures/neel_reversal.svg')
fig.savefig(ROOT/'figures/neel_reversal.png',dpi=180)
print(ROOT/'figures/neel_reversal.svg')
