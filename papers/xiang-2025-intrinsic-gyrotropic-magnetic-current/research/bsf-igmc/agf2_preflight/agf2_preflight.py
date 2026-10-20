"""AgF2 BSF x IGMC: analytic symmetry/DFT input gates (no DFT runs).

Axes are the article's Cartesian x=[100], y=[010], z=[001].
Chi maps a polar current j_i <- axial AC magnetic field B_j. For real
phase-resolved decomposition take eta_C=+1, eta_D=-1, with temporal parity
for an *inter-domain* connector; this does not assert a general fixed-frequency
complex-susceptibility identity. Unit transformations must be consistent.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Any
import hashlib
import json

R_MX = ((-1,0,0),(0,1,0),(0,0,1))
R_MZ = ((1,0,0),(0,1,0),(0,0,-1))
AXES = 'xyz'
CHANNELS = ('C','D')

def _det3(r:tuple) -> int:
    return (r[0][0]*(r[1][1]*r[2][2]-r[1][2]*r[2][1])
            -r[0][1]*(r[1][0]*r[2][2]-r[1][2]*r[2][0])
            +r[0][2]*(r[1][0]*r[2][1]-r[1][1]*r[2][0]))

def response_sign(r:tuple,primed:bool,channel:str,i:int,j:int)->int:
    """Signed diagonal R action on a homogeneous polar<-axial tensor."""
    if channel not in CHANNELS or i not in range(3) or j not in range(3):
        raise ValueError('invalid channel/component')
    if any(r[a][b] for a in range(3) for b in range(3) if a!=b):
        raise ValueError('only diagonal operations supported in this preflight')
    return (_det3(r) * (-1 if channel=='D' and primed else 1)
            *r[i][i]*r[j][j])

def allowed_components(channel:str)->list[str]:
    return [f'{AXES[i]}{AXES[j]}' for i in range(3) for j in range(3)
            if response_sign(R_MX,True,channel,i,j)==+1]

def partner_sign(connector:str,channel:str,component:str)->int:
    if len(component)!=2 or any(a not in AXES for a in component):
        raise ValueError('invalid component')
    if connector=='m_z': return response_sign(R_MZ,False,channel,AXES.index(component[0]),AXES.index(component[1]))
    if connector=="m_z'":return response_sign(R_MZ,True,channel,AXES.index(component[0]),AXES.index(component[1]))
    raise ValueError('unknown connector')

def ahe_z_sign(r:tuple,primed:bool)->int:
    """Sigma_xy^AHE behaves as T-odd axial z under inter-domain mapping."""
    return _det3(r)*r[2][2]*(-1 if primed else 1)

def symmetry_result()->dict[str,Any]:
    return {
      'protocol':'real-phase spin-Zeeman IGMC; polar current <- axial AC B',
      'axes':'x=[100], y=[010], z=[001], Cartesian',
      'source_magnetic_point_group':"m'_x (published for AgF2 Type-I, x-axis SOC easy axis)",
      'single_state_type_I':{
        ch:{'allowed_components':allowed_components(ch),'dimension':len(allowed_components(ch)),
            'xx_allowed':'xx' in allowed_components(ch),'xy_allowed':'xy' in allowed_components(ch)}
        for ch in CHANNELS},
      'switches':{
       kind:{'connector':conn,'C_xy':partner_sign(conn,'C','xy'),
             'D_xx':partner_sign(conn,'D','xx'),
             'AHE_xy':ahe_z_sign(R_MZ,conn.endswith("'")),
             'P_z':R_MZ[2][2]}
       for kind,conn in [('type_I','m_z'),('type_II_strained',"m_z'")]},
      'caveats':[
        'Type-II single-state magnetic group and stable magnetization must be verified from SOC magnetic structure; parity uses literature connector only.',
        'Allowed tensor element may vanish accidentally; conduction IGMC at T=0 vanishes for an insulator.',
        'No Wannier spin/position operator provided; no quantitative real-material IGMC.',
        'C/D antiunitary labels refer to declared real-phase model convention; complex causal response needs independent check.'
      ]}

REQUIRED_PHYSICS=('functional','vdw','hub_u_ag_ev','soc','vacuum_axis','magnetic_axis','magnetic_order','strain_percent','kmesh','ecut_ev','pseudopotential_lock_sha256')
REQUIRED_STRUCTURES=('B1','B2')

def validate_poscar_minimal(path:Path)->list[str]:
    """Validate VASP-5 species/stoichiometry/shape only, not physical provenance."""
    issues=[]
    lines=[x.strip() for x in path.read_text(encoding='utf-8').splitlines() if x.strip()]
    if len(lines)<9:return ['POSCAR_TOO_SHORT']
    try:
        scale=float(lines[1]); assert scale>0
        vec=[[float(x) for x in line.split()] for line in lines[2:5]]
        assert all(len(v)==3 for v in vec)
        species=lines[5].split(); counts=[int(x) for x in lines[6].split()]
        if len(species)!=len(counts) or len(species)!=2 or set(species)!={'Ag','F'} or any(x<=0 for x in counts):
            issues.append('POSCAR_WRONG_SPECIES_OR_COUNTS')
        else:
            atom_count=dict(zip(species,counts))
            if atom_count['F']!=2*atom_count['Ag']:
                issues.append('POSCAR_WRONG_STOICHIOMETRY')
        pos=7
        if lines[pos].lower().startswith('s'):pos+=1
        if not lines[pos].lower().startswith(('d','c','k')):
            issues.append('POSCAR_COORDINATE_MODE_UNRECOGNIZED')
        pos+=1
        if len(lines[pos:])<sum(counts):issues.append('POSCAR_COORDINATES_INCOMPLETE')
        else:
            for row in lines[pos:pos+sum(counts)]:
                parts=row.split()
                if len(parts)<3:
                    issues.append('POSCAR_COORDINATE_TRUNCATED');break
                tuple(float(a) for a in parts[:3])
    except (ValueError,IndexError,AssertionError,KeyError):
        issues.append('POSCAR_PARSE_FAILED')
    return issues


def validate_dft_contract(payload:dict[str,Any],root:Path)->dict[str,Any]:
    """Fail closed before computing: no guessed atoms, U, SOC, POTCAR."""
    reasons=[]
    if payload.get('material')!='AgF2' or payload.get('task')!='bsf_igmc_stage4':
        reasons.append('WRONG_MATERIAL_OR_TASK')
    phi=payload.get('physics') or {}
    for name in REQUIRED_PHYSICS:
        if phi.get(name) in (None,'','UNSET'):
            reasons.append('MISSING_PHYSICS_'+name)
    if phi.get('soc') is not True:reasons.append('SOC_REQUIRED_FOR_IGMC')
    if phi.get('magnetic_axis') not in ('x','y','z'):reasons.append('MAGNETIC_AXIS_UNRESOLVED')
    if phi.get('magnetic_order') not in ('type_I','type_II'):reasons.append('MAGNETIC_CONFIGURATION_UNRESOLVED')
    if phi.get('vacuum_axis')!='z':reasons.append('LAB_AXIS_SETTING_NOT_FROZEN')
    if phi.get('magnetic_order')=='type_I' and phi.get('strain_percent')!=0:
        reasons.append('TYPE_I_SOURCE_STRAIN_SHOULD_BE_ZERO')
    if phi.get('magnetic_order')=='type_II' and phi.get('strain_percent')!=4:
        reasons.append('TYPE_II_SOURCE_STRAIN_SHOULD_BE_FOUR_PERCENT')
    sh=phi.get('pseudopotential_lock_sha256')
    if sh not in (None,'','UNSET') and (not isinstance(sh,str) or len(sh)!=64 or any(c not in '0123456789abcdef' for c in sh.lower())):
        reasons.append('INVALID_POTCAR_HASH')
    structs=payload.get('structures') or {}
    hashes={}
    for key in REQUIRED_STRUCTURES:
        d=structs.get(key) or {}
        name=d.get('poscar')
        digest=d.get('sha256')
        if not name or not digest:
            reasons.append('MISSING_ATOMIC_STRUCTURE_'+key);continue
        if Path(str(name)).is_absolute() or '..' in Path(str(name)).parts:
            reasons.append('INVALID_STRUCTURE_PATH_'+key);continue
        p=(root/name)
        if not p.is_file():reasons.append('STRUCTURE_NOT_FOUND_'+key);continue
        actual=hashlib.sha256(p.read_bytes()).hexdigest()
        if actual!=digest:reasons.append('STRUCTURE_SHA_MISMATCH_'+key)
        else:
            errors=validate_poscar_minimal(p)
            if errors:
                reasons.extend(key+'_'+err for err in errors)
            else: hashes[key]=actual
    if len(hashes)==2 and hashes['B1']==hashes['B2']:
        reasons.append('STRUCTURAL_PARTNERS_IDENTICAL')
    return {'status':'READY_FOR_EXPERT_DFT_INPUT_AUDIT' if not reasons else 'BLOCKED_INPUTS',
            'reasons':reasons, 'verified_structure_sha256':hashes,
            'note':'This is an INPUT-PRECHECK only; passing does not verify physical geometry, magnetic order, relaxation or DFT convergence.'}

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--contract',type=Path,default=Path(__file__).parent/'inputs/agf2_dft_contract.example.json')
    ap.add_argument('--out',type=Path)
    args=ap.parse_args()
    result={'symmetry':symmetry_result(),
            'dft_input_gate':validate_dft_contract(json.loads(args.contract.read_text()),args.contract.parent)}
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(result,indent=2,ensure_ascii=False))
