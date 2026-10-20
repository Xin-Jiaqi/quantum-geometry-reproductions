"""Constrain the microscopic sliding response by m_z and C3.

Results: gapped 8band synthetic KM, equilibrium/static, not driven DC.
"""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
from dataclasses import replace
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from bsf_nmee.sliding_tight_binding import SlidingBilayer, lattice_geometry
from bsf_nmee.symmetry_scaling import split_parity
from bsf_nmee.valley_quadrature import valley_quadrature
from scripts.run_stage3_sign_reversal import spec_at


def sym_irreps(mat):
    tr=np.trace(mat)/2
    nem=(mat+mat.T)/2-tr*np.eye(2)
    anti=(mat-mat.T)/2
    return {'trace':float(2*tr),'isotropic':float(tr),
            'nematic_norm':float(np.linalg.norm(nem)),
            'antisymmetric_norm':float(np.linalg.norm(anti)),
            'matrix':mat.tolist()}


def run(path,quick=False):
    reciprocal=SlidingBilayer().model().reciprocal
    coarse=(29,1.5,25,48) if quick else (41,1.5,35,64)
    fine=(41,1.5,35,64) if quick else (61,1.5,60,80)
    k,w,meta=valley_quadrature(reciprocal,*coarse)
    t=.005
    us=[.025,.05,.1,.2,.4]
    data={}
    for u in us:
        s=replace(spec_at(u),tunneling_eV=t)
        result=split_parity(s,k,w)
        data[f'{u:.3f}']={key:sym_irreps(result[key]) for key in ('full','even','mixed_odd')}
        data[f'{u:.3f}']['min_mu_distance_eV']=result['mu_distance']
        print('u',u,'trace/u³',data[f'{u:.3f}']['full']['trace']/u**3,
              'mixed_nematic/u',data[f'{u:.3f}']['mixed_odd']['nematic_norm']/u,flush=True)
    a,_=lattice_geometry()
    rotations={}
    for angle in (60,120,180):
        rad=np.deg2rad(angle); r=np.array([[np.cos(rad),-np.sin(rad)],[np.sin(rad),np.cos(rad)]])
        s=replace(spec_at(.1),tunneling_eV=t)
        shift=(r@(np.array(s.shift)@a))@np.linalg.inv(a)
        calc=split_parity(replace(s,shift=tuple(shift)),k,w)
        base=np.asarray(data['0.100']['full']['matrix'])
        if angle==120:expected=r@base@r.T
        elif angle==180:expected=-base
        else:expected=None
        rotations[str(angle)]={'fractional_shift':shift.tolist(),
          'full':sym_irreps(calc['full']),
          'predicted_full_matrix':None if expected is None else expected.tolist(),
          'matrix_residual':None if expected is None else float(np.max(abs(calc['full']-expected))),
          'trace_expected_ratio':-1 if angle in (60,180) else 1,
          'trace_ratio':float(np.trace(calc['full'])/np.trace(base))}
    comparisons={}
    k2,w2,meta2=valley_quadrature(reciprocal,*fine)
    for u in (.05,.1,.2):
        s=replace(spec_at(u),tunneling_eV=t)
        c=split_parity(s,k,w)
        f=split_parity(s,k2,w2)
        b=split_parity(s,k,w,step_T=.125)
        comparisons[f'{u:.3f}']={'coarse_trace':float(np.trace(c['full'])),
          'fine_trace':float(np.trace(f['full'])),
          'fine_minus_coarse_trace':float(np.trace(f['full']-c['full'])),
          'max_component_quadrature_delta':float(np.max(abs(f['full']-c['full']))),
          'max_component_field_delta':float(np.max(abs(b['full']-c['full']))),
          'fine_mixed_nematic_norm':sym_irreps(f['mixed_odd'])['nematic_norm']}
        print('mesh check',u, comparisons[f'{u:.3f}'],flush=True)
    # Scale ratios as empirical evidence, not symmetry proof of individual channels.
    trace_ratio_2=float(data['0.100']['full']['trace']/data['0.050']['full']['trace'])
    nem_ratio_2=float(data['0.100']['mixed_odd']['nematic_norm']/data['0.050']['mixed_odd']['nematic_norm'])
    blob={'scope':'symmetry of 8-band microscopic sliding, 0 K, static equilibrium response',
        'limits':'no DFT/Wannier material; no validated driven DC; no true ferroelectric switching energy',
        'model':'C3-invariant unshifted homobilayer and layer-exchanging horizontal mirror m_z',
        'indices':'matrix spin x/y; B_z; electric x/y',
        'tunneling_eV':t,'coarse':{'grid':coarse,'meta':meta},
        'fine':{'grid':fine,'meta':meta2},
        'data':data,'rotated_shifts':rotations,'convergence':comparisons,
        'empirical_scale_ratios_u01_over_u005':{'trace':trace_ratio_2,'mixed_nematic':nem_ratio_2},
        'theorem':'m_z: M(-tau)=-M(tau); C3: M(Rtau)=R M(tau) R^T; scalar trace has no constant or linear term and begins at cubic if analytic; traceless sector can be linear'}
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(blob,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    return blob

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--quick',action='store_true');p.add_argument('--output',type=Path,default=ROOT/'results/sliding_stage6_symmetry.json');args=p.parse_args()
    run(args.output,quick=args.quick)
