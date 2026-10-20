"""Magnetic-pair IGMC screening from pinned group-theory toolkit generators.

Used as a dimension/selection-rule screen, NOT a complete crystallographic
stacking enumerator. Spatial setting rotations are explicit.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np

SNAPSHOT=Path(__file__).resolve().parents[2]/'data'/'toolkit_mpg_generators_v018.txt'
P_XZ=np.array([[0,0,1],[0,1,0],[1,0,0]],float)

def toolkit_generator(group_number:int):
    if group_number not in (2,6,8,9,11):raise ValueError('restricted, audited MPG IDs')
    mapping={}
    records={}
    for line in SNAPSHOT.read_text(encoding='utf-8').splitlines():
        if line.startswith('MAP '):
            family,content=line[4:].split(' ',1)
            mapping[family]={}
            for token in content.split(';'):
                code,name,flat=token.split(':',2)
                mapping[family][code]=np.array([float(x) for x in flat.split(',')]).reshape(3,3)
        if line.startswith('RECORDS '):
            for record in line[8:].split('|'):
                head,codes=record.split(':',1)
                idx,hm=head.split('=',1)
                records[int(idx)]=(hm,codes.split(',') if codes else [])
    hm,codes=records[group_number]
    assert len(codes)==1
    key=codes[0];base=key.rstrip('!')
    family=('cubic_Oh' if group_number<=5 else 'tetragonal_D4h')
    R=mapping[family][base]
    return hm,R,key.endswith('!')

def get_op(name):
    ids={'1prime':2,'2z':6,'2zprime':8,'mz':9,'mzprime':11}
    if name in ('mx','mxprime','2x','2xprime'):
        alias={'mx':'mz','mxprime':'mzprime','2x':'2z','2xprime':'2zprime'}[name]
        hm,R,tr=toolkit_generator(ids[alias]);return P_XZ@R@P_XZ.T,tr
    if name=='identity':return np.eye(3),False
    hm,R,tr=toolkit_generator(ids[name]);return R,tr

def tensor_transform(R,tr,channel):
    if channel not in ('C','D'):raise ValueError(channel)
    eta=1 if channel=='C' else -1
    sign=(eta if tr else 1)*round(np.linalg.det(R))
    return np.column_stack([(sign*R@b.reshape(3,3)@R.T).ravel() for b in np.eye(9)])

def basis_for_single_state(name,channel):
    if name=='identity':return np.eye(9)
    R,tr=get_op(name)
    _,sing,vh=np.linalg.svd(np.eye(9)-tensor_transform(R,tr,channel),full_matrices=True)
    return vh[np.sum(sing>1.e-8):,:].T

def allowed(basis,i,j):
    return bool(basis.shape[1]>0 and np.max(np.abs(basis[3*i+j,:]))>1e-8)

def switching_sign(connector,channel,i,j):
    R,tr=get_op(connector)
    unit=np.zeros((3,3));unit[i,j]=1
    out=((1 if channel=='C' else -1) if tr else 1)*round(np.linalg.det(R))*R@unit@R.T
    if np.max(np.abs(out-np.sign(out[i,j])*unit))>1e-10:
        return None # component mixes; cannot assert scalar switch parity
    return int(np.sign(out[i,j]))

def screen():
    groups={}
    for name in ('identity','1prime','mx','mxprime'):
        groups[name]={}
        for channel in ('C','D'):
            b=basis_for_single_state(name,channel)
            groups[name][channel]={'dimension':int(b.shape[1]),
                                   'xy_allowed':allowed(b,0,1),
                                   'xx_allowed':allowed(b,0,0)}
    pairs={}
    for name in ('mz','mzprime','2x','2xprime'):
        R,tr=get_op(name)
        pairs[name]={'time_reversal':tr,'polar_Pz_sign':int(round(R[2,2])),
                     'Cxy_switch_sign':switching_sign(name,'C',0,1),
                     'Dxx_switch_sign':switching_sign(name,'D',0,0)}
    return {'source':'Xin-Jiaqi/group-theory-operations-toolkit v0.18.0 pinned generators',
            'setting':'Cartesian x,y,z; mx = P_xz mz P_xz^T',
            'internal_state_groups':groups,'connection_pairs':pairs}
