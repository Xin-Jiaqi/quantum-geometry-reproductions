"""Full AB→BA path, local polar-stack expansion and microscopic null controls."""
import json
from dataclasses import asdict,replace
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from bsf_nmee.bilayer_audit import AuditSpec,AuditBilayer,integration_grid,symmetry_residuals
from bsf_nmee.bilayer_protocol import ab_ba_shift,half_filled_spec
from bsf_nmee.static_observables import observables


def main():
    k,w,meta=integration_grid(AuditBilayer(),61,60,96)
    result={'scope':'half-filled gapped mirror-paired model; no structural energetics',
            'grid':meta,'coordinates':'60deg code; A1=a1,A2=a2-a1 is the 120deg convention'}
    def evaluate(spec):
        try:
            spec=half_filled_spec(spec,k)
            return {'parameters':asdict(spec),**observables(AuditBilayer(spec),k,w,step_T=.0625)}
        except ValueError as e:
            return {'parameters':asdict(spec),'excluded':str(e)}
    def save():
        (ROOT/'results/bilayer_stage9_path.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    result['path']=[]
    for t in (.005,.02):
        for lam in np.linspace(0,1,17):
            row={'lambda':float(lam),'t_eV':t,**evaluate(AuditSpec(shift=ab_ba_shift(lam),tunneling_eV=t))}
            result['path'].append(row)
            print('path',t,lam,row.get('trace'),row.get('difference'),flush=True)
        save()
    result['polar_local']=[]
    ab=np.array(ab_ba_shift(0))
    for direction in ((.23,.41),(-.41,.23)):
        for q in (0.,-.04,-.02,-.01,.01,.02,.04):
            result['polar_local'].append({'direction':direction,'q':q,
                **evaluate(AuditSpec(shift=tuple(ab+q*np.array(direction)),tunneling_eV=.005))})
            print('local',direction,q,flush=True)
    save()
    result['essential_terms']={}
    reference=AuditSpec(shift=(.023,.041),tunneling_eV=.005)
    for name,changes in {
        'reference':{},'no_Rashba':{'rashba_soc_eV':0.},
        'Rashba_reversed':{'rashba_soc_eV':-.02},'no_staggered':{'staggered_eV':0.},
        'staggered_reversed':{'staggered_eV':-.01},'no_intrinsic_SOC':{'intrinsic_soc_eV':0.},
        'same_only':{'cross_scale':0.},'cross_only':{'same_scale':0.},
        'decoupled':{'same_scale':0.,'cross_scale':0.}}.items():
        spec=replace(reference,**changes)
        result['essential_terms'][name]=evaluate(spec)
        result['essential_terms'][name]['symmetry']=symmetry_residuals(spec)
        print('essential',name,flush=True)
    save()
    # Control the local expansions at AB and both full-path endpoints with finer integration.
    kf,wf,_=integration_grid(AuditBilayer(),101,100,144)
    result['fine']={}
    for name,shift in [('AB',ab_ba_shift(0)),('mid',ab_ba_shift(.5)),('BA',ab_ba_shift(1)),
                       ('AB_plus',tuple(ab+.01*np.array([.23,.41]))),
                       ('AB_minus',tuple(ab-.01*np.array([.23,.41])))]:
        spec=half_filled_spec(AuditSpec(shift=shift,tunneling_eV=.005),kf)
        result['fine'][name]={'parameters':asdict(spec),**observables(AuditBilayer(spec),kf,wf,step_T=.0625)}
        print('fine',name,flush=True)
    save()


if __name__=='__main__':main()
