"""Render plots from the frozen stage-6 JSON without rerunning the model."""
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
r=json.loads((ROOT/'results/sliding_stage6_symmetry.json').read_text())
us=sorted(float(x) for x in r['data'])
full=[r['data'][f'{u:.3f}']['full'] for u in us]
mix=[r['data'][f'{u:.3f}']['mixed_odd'] for u in us]
trace=np.array([item['trace'] for item in full])
nem=np.array([item['nematic_norm'] for item in mix])
fig,ax=plt.subplots(figsize=(6.5,4.2))
ax.plot(np.array(us)**3,trace,'o-',label=r'$\mathrm{Tr}\,\Gamma_\parallel$')
c=trace[0]/us[0]**3
ax.plot(np.array(us)**3,c*np.array(us)**3,'--',label=r'Leading $u^3$ (near origin)')
ax.axhline(0,lw=.5)
ax.set_xlabel(r'Cubed sliding coordinate $u^3$')
ax.set_ylabel(r'$\Gamma_{xz;x}+\Gamma_{yz;y}$ (model units)')
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(ROOT/'results/sliding_stage6_trace_cubic.png',dpi=190)
fig.savefig(ROOT/'results/sliding_stage6_trace_cubic.svg')
plt.close(fig)
fig,ax=plt.subplots(figsize=(6.5,4.2))
ax.plot(us,nem,'o-',label='Mixed-channel nematic response')
ax.plot(us,(nem[0]/us[0])*np.array(us),'--',label=r'Leading $u$ (near origin)')
ax.set_xlabel('Sliding coordinate u')
ax.set_ylabel('Nematic norm (model units)')
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(ROOT/'results/sliding_stage6_nematic_linear.png',dpi=190)
fig.savefig(ROOT/'results/sliding_stage6_nematic_linear.svg')
plt.close(fig)
