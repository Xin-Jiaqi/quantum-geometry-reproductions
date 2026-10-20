#!/usr/bin/env python3
"""Stdlib-only independent Stage-2 C2DB-derived geometry audit. No DFT runs."""
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
STAGE1 = ROOT.parent / 'stage6_harness' / 'stage1_review' / 'structures'
SOURCE_SHA = {
  'B1':'b3c1b942df69b5442e52d0d38fe839b4b61feff3619b2adb4736faa82501c122',
  'B2':'ab2f3d84af0a074079824b1f0b94d928d16609e889ac7b1ba3fd9a50e2c69087',
}
GAPS=(2.6,3.0,3.4)

def load(path):
    l=path.read_text().splitlines()
    assert abs(float(l[1])-1.)<1e-12 and l[5].split()==['Ag','F'] and l[6].split()==['4','8']
    assert l[7].lower().startswith('cart')
    cell=[[float(x) for x in l[i].split()] for i in (2,3,4)]
    for i in range(3):
        assert all(abs(cell[i][j])<1e-10 for j in range(3) if j!=i)
    return [cell[i][i] for i in range(3)],[[float(x) for x in l[i].split()] for i in range(8,20)]

def dist(a,b):
    return abs(((a-b+1.5)%1.)-.5)

def verify():
    m=json.loads((ROOT/'candidates'/'GEOMETRY_MANIFEST.json').read_text())
    assert m['data_origin']=='NOT_RUN' and m['physics_approved_by_user'] is False
    assert m['potcar_lock_sha256'] is None and m['dft_job_submitted'] is False
    assert len(m['cases'])==6
    original={}
    for label in ('B1','B2'):
        p=STAGE1/f'constructed_{label}_candidate'/'POSCAR'
        assert hashlib.sha256(p.read_bytes()).hexdigest()==SOURCE_SHA[label]
        original[label]=load(p)
    maximum=0.
    for gap in GAPS:
        trial={}
        for label in ('B1','B2'):
            path=ROOT/'candidates'/f'{label}_g{gap:.1f}'/'POSCAR'
            lengths,pts=load(path)
            old,oldpts=original[label]
            d=gap-3.0
            assert max(abs(lengths[i]-(old[i]+(d if i==2 else 0))) for i in range(3))<1e-7
            assert len(pts)==12
            for (x,y,z),(X,Y,Z) in zip(oldpts,pts):
                assert max(abs(x-X),abs(y-Y),abs(z+(d if z>old[2]/2 else 0)-Z))<1e-7
            aa,bb,cc=lengths
            assert abs(cc-(max(p[2] for p in pts)-min(p[2] for p in pts))-20.)<1e-6
            low=[p[2] for i,p in enumerate(pts) if i>=4 and p[2]<cc/2]
            high=[p[2] for i,p in enumerate(pts) if i>=4 and p[2]>cc/2]
            assert abs(min(high)-max(low)-gap)<1e-6
            trial[label]=(lengths,pts)
        cell,bb=trial['B2'];a,b,c=cell
        aa=trial['B1'][1]
        used=set()
        for i,(x,y,z) in enumerate(aa):
            target=[(x/a+.5)%1,(y/b+.16)%1,(1-z/c)%1]
            opt=[]
            for j,(X,Y,Z) in enumerate(bb):
                if j in used or (i<4)!=(j<4): continue
                err=max(dist(u,v) for u,v in zip(target,[X/a,Y/b,Z/c]))
                opt.append((err,j))
            assert opt
            err,j=min(opt)
            assert err<3e-6
            maximum=max(maximum,err);used.add(j)
        assert len(used)==12
    print('Stage2 GEOMETRY-ONLY PASS: 6 POSCARs, 20 A vacuum, partner bijection; max fractional error %.2e'%maximum)
    print('DFT=NOT_RUN; POTCAR=BLOCKED; physics approval=FALSE')
    return True

if __name__=='__main__':
    verify()
