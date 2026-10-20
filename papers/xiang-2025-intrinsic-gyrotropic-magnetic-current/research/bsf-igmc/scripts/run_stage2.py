"""Generate model-level evidence without changing the paper reproduction."""
import argparse, json
from pathlib import Path
from dataclasses import asdict, replace
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from bsf_igmc_micro.model import (Parameters,connection_error,geometric_response,
  layer_polarization,wilson_zeeman_slope,band_energy)
from bsf_igmc_micro.symmetry import screen

def build(selected=None):
  p=Parameters(mass=.85)
  out={'scientific_status':'MODEL_VALIDATED_BSF_NOT_ESTABLISHED',
       'model':'8 bands: 2 layers x 2 displaced orbitals x physical spin',
       'onsite_layer_polar_term':0,
       'hopping_source':'explicit real-space Gaussian of layer/orbital center displacement',
       'parameters':asdict(p),
       'conventions':'lattice momentum a=1, e=hbar=1, Dxx/omega reduced BZ averaged kernel (not SI current)',
       'toolkit':screen(),
       'branches':{}}
  for branch,mu in [('AF',0.),('FM',-.03)]:
    if selected is not None and branch!=selected:continue
    b={'mu_eV':mu,'symmetry_errors':connection_error(.19,branch,p,N=5),
       'plus':{},'minus':{},'u0':{},'wilson':{},'ablation':{},'band_energy_only':{}}
    for N in (21,31):
      b['plus'][str(N)]=geometric_response(.19,branch,p,N=N,mu=mu)
    b['minus']['31']=geometric_response(-.19,branch,p,N=31,mu=mu)
    b['u0']['21']=geometric_response(0.,branch,p,N=21,mu=mu)
    for u in [0.,.19,-.19]:
      b[f'layer_charge_{u:+.2f}']=layer_polarization(u,branch,p,N=31,mu=mu)
    for N in (21,):
      b['wilson'][str(N)]=wilson_zeeman_slope(.19,branch,p,N=N,field=3e-4)
    for name,val in [('inter_spin',0.),('inter_soc',0.),('lambda_so',0.),('cant',0.)]:
      pb=replace(p,**{name:val})
      b['ablation'][name]=geometric_response(.19,branch,pb,N=21,mu=mu)
    for u in (0.,.19,.38,.48):
      b['band_energy_only'][f'{u:.3f}']=band_energy(u,branch,p,N=21)
    # A structural double well is NOT validated without ionic & electrostatic energy.
    out['branches'][branch]=b
  return out

if __name__=='__main__':
  parser=argparse.ArgumentParser()
  parser.add_argument('--output',default='results/stage2_eightband.json')
  parser.add_argument('--branch',choices=['AF','FM'])
  args=parser.parse_args()
  v=build(args.branch)
  target=Path(args.output);target.parent.mkdir(exist_ok=True,parents=True)
  target.write_text(json.dumps(v,indent=2,ensure_ascii=False)+'\n')
  for name,data in v['branches'].items():
    print(name,'D(+u)=',data['plus']['31']['Dxx_over_omega'],
        'D(-u)=',data['minus']['31']['Dxx_over_omega'],
        'gap>=',data['plus']['31']['indirect_gap'],
        'Pz(+u)=',data['layer_charge_+0.19'],
        'Wilson(31)=',data['wilson']['21'])
  print('Output:',target)
