"""Figures from frozen stage-5 data; no calculations are run by plotting."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[1]
d=json.loads((root/'results/sliding_stage5_weak_coupling.json').read_text())
fig,ax=plt.subplots(figsize=(7.8,4.7))
for u,content in d['data'].items():
    rr=content['sweep']
    ax.plot([r['t_eV']*1e3 for r in rr],[r['plus']*1e6 for r in rr],marker='o',label=f'slide u={u}')
ax.axhline(0,color='gray',lw=.8)
ax.set_xlabel('Interlayer hopping scale (meV)')
ax.set_ylabel(r'Total static $\Gamma_{xz;x}$ ($10^{-6}$ in code units)')
ax.set_title('A sign reversal at fixed sliding, without a sampled gap closure')
ax.grid(alpha=.18);ax.legend()
fig.tight_layout()
fig.savefig(root/'results/sliding_stage5_hopping_transition.png',dpi=190)
fig.savefig(root/'results/sliding_stage5_hopping_transition.svg');plt.close(fig)
fig,ax=plt.subplots(figsize=(7.8,4.7))
us=list(d['data']); even=[d['data'][u]['sweep'][-1]['even_under_cross_sign']*1e6 for u in us]
odd=[d['data'][u]['sweep'][-1]['odd_under_cross_sign']*1e6 for u in us]
ax.plot([float(u) for u in us],even,'o-',label='Even under cross-hopping sign')
ax.plot([float(u) for u in us],odd,'s-',label='Odd under cross-hopping sign')
ax.plot([float(u) for u in us],[even[i]+odd[i] for i in range(len(us))],'^-',label='Total')
ax.axhline(0,color='gray',lw=.8)
ax.set_xlabel('Sliding coordinate u');ax.set_ylabel(r'Static $\Gamma_{xz;x}$ ($10^{-6}$ in code units)')
ax.set_title('Mechanism: positive even contribution overtakes negative interference')
ax.grid(alpha=.18);ax.legend()
fig.tight_layout();fig.savefig(root/'results/sliding_stage5_parity_competition.png',dpi=190)
fig.savefig(root/'results/sliding_stage5_parity_competition.svg');plt.close(fig)
