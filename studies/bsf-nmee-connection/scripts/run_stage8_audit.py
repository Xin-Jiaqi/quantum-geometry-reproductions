"""Reproducible local physics audit; run from repo root with declared PYTHONPATH.

No DFT/material prediction. All cases use the same static gapped EB protocol.
"""
from dataclasses import asdict, replace
from datetime import datetime, timezone
from pathlib import Path
import argparse
import json
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from bsf_nmee.bilayer_audit import AuditSpec, AuditBilayer, integration_grid, symmetry_residuals
from bsf_nmee.static_observables import observables


def at(u, angle=0., **changes):
    spec = AuditSpec(shift=(.23*u,.41*u),tunneling_eV=.005,**changes)
    if angle:
        model=AuditBilayer(spec); theta=np.deg2rad(angle)
        r=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
        spec=replace(spec,shift=tuple(np.array(spec.shift)@model.a@r.T@np.linalg.inv(model.a)))
    return spec


def run(output, quick=False):
    coarse=(29,25,48) if quick else (61,60,96)
    fine=(41,35,64) if quick else (101,100,144)
    k,w,meta=integration_grid(AuditBilayer(),*coarse)
    kf,wf,metaf=integration_grid(AuditBilayer(),*fine)
    def evaluate(spec, *, grid=None, step=.125, flat=None):
        kk,ww=(k,w) if grid is None else grid
        result=observables(AuditBilayer(spec),kk,ww,step_T=step,flatten_gap_eV=flat)
        return {'parameters':asdict(spec),'step_T':step,'flatten_gap_eV':flat,**result}
    result={
        'scope':'synthetic 8-band bilayer; gapped T=0; static d_Bz spin-electric susceptibility',
        'created_utc':datetime.now(timezone.utc).isoformat(),
        'grid':{'coarse':coarse,'fine':fine,'coarse_meta':meta,'fine_meta':metaf,
                'C3_closed':True,'normalization':'d^2k/(2pi)^2; no fitted or projected observables'},
        'units':{'C':'legacy effective-g spin-density coefficient, 1/(nm V T)',
                 'dipole_over_e_d_per_cell':'D_z/(positive e * declared layer separation); dimensionless',
                 'dipole_over_e_d_per_nm2':'P_z/(positive e * layer separation), nm^-2',
                 'u':'dimensionless; fractional displacement u*(0.23,0.41)'},
        'limitations':['not driven DC','not a material Wannier model','no ionic relaxation or energy barrier',
                       'channel controls change Hamiltonians; not unique contributions to one observable',
                       'flattening is an artificial equal-spectrum comparison, not an experimental knob'],
    }
    def save():
        output.parent.mkdir(parents=True,exist_ok=True)
        output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    result['symmetry']={name:symmetry_residuals(spec) for name,spec in {
        'mirror_default':at(.1),
        'C3_broken':at(.1,bond_anisotropy_eV=.008),
        'same_Rashba_direct':at(.1,stacking='direct'),
        'graphene_direct':at(.1,stacking='direct',staggered_eV=0.,rashba_soc_eV=0.,mu_eV=0.)}.items()}
    result['sliding']=[]
    for u in (0.,.025,.05,.1,.2,.4,.8,1.):
        point={'u':u,**evaluate(at(u))};result['sliding'].append(point)
        print('slide',u,'dipole',point['dipole_over_e_d_per_cell'],
              'trace',point['trace'],'diff',point['difference'],flush=True)
    save()
    result['convergence']=[]
    for u in (.025,.1,.4):
        result['convergence'].append({'u':u,'fine':evaluate(at(u),grid=(kf,wf)),
                                     'half_B':evaluate(at(u),step=.0625)})
        print('convergence',u,flush=True)
    result['angles']=[]
    for angle in (0.,30.,60.,90.,120.,180.):
        row={'angle_deg':angle,'small':evaluate(at(.05,angle)),
             'large':evaluate(at(.1,angle))}
        result['angles'].append(row)
        print('angle',angle,flush=True)
    save()
    result['parameter_controls']={}
    variations={'rashba_minus20':{'rashba_soc_eV':.016},
                'rashba_plus20':{'rashba_soc_eV':.024},
                'intrinsic_minus20':{'intrinsic_soc_eV':.008},
                'intrinsic_plus20':{'intrinsic_soc_eV':.012},
                'staggered_minus20':{'staggered_eV':.008},
                'staggered_plus20':{'staggered_eV':.012},
                'C3_broken':{'bond_anisotropy_eV':.008},
                'same_Rashba_direct':{'stacking':'direct'},
                'graphene_direct':{'stacking':'direct','staggered_eV':0.,'rashba_soc_eV':0.,'mu_eV':0.}}
    for name,changes in variations.items():
        rows=[]
        for u in (.025,.05,.1):
            try:
                rows.append({'u':u,**evaluate(at(u,**changes))})
            except ValueError as error:
                rows.append({'u':u,'excluded':str(error),'parameters':asdict(at(u,**changes))})
        result['parameter_controls'][name]=rows
        print('parameters',name,flush=True)
    save()
    result['channels']=[]
    for t in (.001,.002,.005,.01):
        spec=replace(at(.2),tunneling_eV=t)
        row={'t_eV':t}
        for name,s,c in [('full',1,1),('same_only',1,0),('cross_only',0,1),
                         ('relative_sign',1,-1),('decoupled',0,0)]:
            row[name]=evaluate(replace(spec,same_scale=s,cross_scale=c))
        result['channels'].append(row)
        print('channels',t,flush=True)
    result['flat_spectrum']=[]
    for t in (.005,.02):
        for u in (.025,.05,.1,.2,.4):
            spec=replace(at(u),tunneling_eV=t)
            result['flat_spectrum'].append({'u':u,'t_eV':t,'physical':evaluate(spec),
                                            'flat':evaluate(spec,flat=.03)})
            print('flat',t,u,flush=True)
    save()
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--quick',action='store_true')
    parser.add_argument('--output',type=Path,default=ROOT/'results/bilayer_stage8_audit.json')
    args=parser.parse_args()
    run(args.output,args.quick)
