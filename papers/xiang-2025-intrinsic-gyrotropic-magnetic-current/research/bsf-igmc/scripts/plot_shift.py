"""Figures from independently calculated toy-model scan (optional matplotlib)."""
from pathlib import Path
import csv
import matplotlib.pyplot as plt
src=Path('results/shift_scan_17.csv')
rows=list(csv.DictReader(src.open()))
for value,label,filename in (
  ('Dxx_over_omega',r'IGMC geometric kernel $D_{xx}/\omega$ (reduced units)','results/shift_Dxx_8band.png'),
  ('layer_charge',r'Layer charge imbalance $\langle\tau_z\rangle$ (electrons/cell)','results/shift_layer_charge_8band.png')):
  fig,ax=plt.subplots(figsize=(6.8,4.4))
  for branch,display in [('AF',r'unitary $m_z$'),('FM',r'antiunitary $m_z\mathcal{T}$')]:
    sub=[r for r in rows if r['branch']==branch]
    ax.plot([float(r['u']) for r in sub],[float(r[value]) for r in sub],marker='o',markersize=3,label=display)
  ax.axhline(0,lw=.8,color='black',alpha=.5)
  ax.set_xlabel('Relative layer displacement u / lattice period')
  ax.set_ylabel(label)
  ax.set_title('8-band bilayer toy model: no imposed layer-odd onsite term')
  ax.grid(alpha=.22);ax.legend();fig.tight_layout()
  fig.savefig(filename,dpi=160);plt.close(fig)
  print(filename)
