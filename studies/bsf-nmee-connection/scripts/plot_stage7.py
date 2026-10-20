"""Plots of frozen stage-7 valley k.p numerical results (not fit to PRL)."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/'results'

capture=json.loads((D/'sliding_stage7_capture.json').read_text())
fig,ax=plt.subplots(figsize=(7,4.7))
for u in (.1,.4):
  for t in (.005,.02):
    items=[v for v in capture if v['u']==u and v['t_eV']==t]
    items.sort(key=lambda r:r['radius_inv_nm'])
    ax.plot([v['radius_inv_nm'] for v in items],
            [100*abs(v['disk_fraction_signed_xx']-1) for v in items],
            marker='o',label=f'u={u}, t={1000*t:g} meV')
ax.set_yscale('log');ax.set_xlabel('Valley disk radius (nm$^{-1}$)')
ax.set_ylabel('Signed-valley vs whole BZ |ratio - 1| (%)')
ax.set_title('K/K\u2032 disks capture nearly all signed static NMEE')
ax.legend(fontsize=8);fig.tight_layout();fig.savefig(D/'sliding_stage7_valley_capture.png',dpi=170);fig.savefig(D/'sliding_stage7_valley_capture.svg');plt.close(fig)

channels=json.loads((D/'sliding_stage7_channels.json').read_text())
fig,ax=plt.subplots(figsize=(7,4.7))
u=np.array([r['u'] for r in channels]);
for label,key,mark in [('BZ mixed','BZ_mixed_coeff','o'),('Valley k.p mixed','disk_kp2_mixed_coeff','x'),
                       ('BZ even','BZ_even_coeff','s'),('Valley k.p even','disk_kp2_even_coeff','+')]:
    ax.plot(u,[r[key] for r in channels],marker=mark,label=label)
ax.axhline(0,color='0.3',linewidth=.8)
ax.set_xlabel('Sliding path parameter u')
ax.set_ylabel('Quadratic response coefficient (program units / eV$^2$)')
ax.set_title('Low-energy valleys reproduce channel competition')
ax.legend(fontsize=9);fig.tight_layout();fig.savefig(D/'sliding_stage7_valley_channels.png',dpi=170);fig.savefig(D/'sliding_stage7_valley_channels.svg');plt.close(fig)

roots=json.loads((D/'sliding_stage7_root_and_trace.json').read_text())
fig,ax=plt.subplots(figsize=(7,4.7))
x=[r['u'] for r in roots['root_scan']]
ax.plot(x,[r['full_disks_xx'] for r in roots['root_scan']],marker='o',label='Full 8-band, K/K\u2032 patch')
ax.plot(x,[r['kp2_disks_xx'] for r in roots['root_scan']],linestyle='--',marker='x',label='Quadratic valley k.p')
ax.axhline(0,color='0.3',linewidth=.8)
ax.axvline(roots['estimated_zero_kp2'],linestyle=':',label='Valley k.p zero')
ax.set_xlabel('Sliding path parameter u')
ax.set_ylabel('Static $\\Gamma_{xz;x}$ (program units)')
ax.set_title('Valley k.p reproduces sliding sign reversal')
ax.legend(fontsize=8);fig.tight_layout();fig.savefig(D/'sliding_stage7_valley_sign_reversal.png',dpi=170);fig.savefig(D/'sliding_stage7_valley_sign_reversal.svg');plt.close(fig)
