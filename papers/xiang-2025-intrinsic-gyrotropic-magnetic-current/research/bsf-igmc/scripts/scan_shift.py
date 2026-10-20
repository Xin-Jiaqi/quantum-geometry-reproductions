"""Exploratory u dependence at fixed parameters; output not claimed as a FE path."""
import csv,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from bsf_igmc_micro.model import Parameters,geometric_response,layer_polarization

def scan(N=17):
  p=Parameters(mass=.85)
  rows=[]
  for b,mu in [('AF',0.),('FM',-.03)]:
    for i in range(-6,7):
      u=i*.04
      r=geometric_response(u,b,p,N=N,mu=mu,strict_gap=False)
      rows.append({'branch':b,'u':round(u,5),'Dxx_over_omega':r['Dxx_over_omega'],
         'layer_charge':layer_polarization(u,b,p,N=N,mu=mu),
         'indirect_gap_eV':r['indirect_gap'], 'occupied_band_range':str(r['occ_range']),
         'mesh':N,'mu_eV':mu})
  return rows

if __name__=='__main__':
  import argparse
  pa=argparse.ArgumentParser();pa.add_argument('--output',default='results/shift_scan_17.csv');pa.add_argument('--mesh',type=int,default=17);a=pa.parse_args()
  rows=scan(N=a.mesh)
  path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
  with path.open('w',newline='') as h:
    w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
  print('saved',path,len(rows),'rows')
