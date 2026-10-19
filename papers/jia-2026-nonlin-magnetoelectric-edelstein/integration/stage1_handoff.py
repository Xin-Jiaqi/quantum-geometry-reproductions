"""Strict, read-only preflight for existing W90Flow Stage-1 / NLO Stage-2 package.

Never assumes Wannier HR+centres are sufficient for quantitative NMEE.
The input SHA256 contract matches nlo-response-workflow manifest.py.
"""
from __future__ import annotations
import hashlib
import json
import sys
from pathlib import Path
REQUIRED=('POSCAR','wannier90_hr.dat','wannier90-wcenter.dat','quality.json')

def _sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1<<20),b''):h.update(block)
    return h.hexdigest()

def _load(path):
    v=json.loads(Path(path).read_text(encoding='utf-8'))
    if not isinstance(v,dict):raise ValueError(f'JSON must be object: {path}')
    return v

def inspect_package(root,*,structure=None):
    root=Path(root).resolve();f=root/'package-manifest.json'
    m=_load(f); ss=m.get('structures')
    if not isinstance(ss,dict) or not ss:raise ValueError('missing structures')
    if structure is None:
        if len(ss)!=1:raise ValueError('multiple structures; choose --structure')
        structure=next(iter(ss))
    if structure not in ss:raise ValueError('unknown structure')
    members=ss[structure]
    if not isinstance(members,dict) or not isinstance(members.get('files'),dict):
        raise ValueError('structure.files must be a mapping')
    files=members['files']; where=root/structure
    hashes={}
    for name in REQUIRED:
        if name not in files or not isinstance(files[name],str) or len(files[name])!=64:
            raise ValueError(f'missing SHA256: {name}')
        path=where/name
        if not path.is_file():raise ValueError(f'missing {path}')
        value=_sha(path)
        if value!=files[name]:raise ValueError(f'SHA256 mismatch: {name}')
        hashes[name]=value
    q=_load(where/'quality.json')
    if q.get('passed') is not True:raise ValueError('Wannier quality gate not passed')
    nw=int(m.get('num_wann',0)); no=int(m.get('nocc_model',0))
    if not 0<no<nw:raise ValueError('invalid nocc_model')
    if int(q.get('num_wann',-1))!=nw or int(q.get('nocc_model',-1))!=no:
        raise ValueError('quality/model dimensions disagree')
    with open(where/'wannier90_hr.dat') as hr:
        hr.readline();nw_hr=int(hr.readline().strip())
    if nw_hr!=nw:raise ValueError('HR dimension mismatch')
    centres=sum(1 for line in (where/'wannier90-wcenter.dat').read_text().splitlines()
                if len(line.split())>=4 and line.split()[0]=='X')
    if centres!=nw:raise ValueError('Wannier centre count mismatch')
    phys=m.get('physics',{}) or {}
    if not isinstance(phys,dict):raise ValueError('physics must be a mapping')
    d=phys.get('dimensionality')
    if d!=2:raise ValueError('NMEE 2D research backend requires explicit dimensionality=2')
    soc=phys.get('soc',q.get('soc'))
    if not isinstance(soc,bool):raise ValueError('explicit SOC bool required for spin basis provenance')
    if isinstance(q.get('soc'),bool) and q['soc']!=soc:
        raise ValueError('manifest and quality SOC conflict')
    return {'handoff_contract':'W90Flow_topK_package_v1','structure':structure,
            'status':'BLOCKED_PHYSICAL_OPERATORS','sha256_verified':True,
            'manifest_sha256':_sha(f),'files_sha256':hashes,'num_wann':nw,'nocc_model':no,
            'soc':soc,'dimensionality':d,
            'missing_required_for_NMEE':['same-gauge projected spin/Zeeman operators S(k)',
            'projected Berry/position connection A_beta(k), or an explicitly restricted approximation',
            'consistent k derivatives for projected spin and symmetry validation'],
            'warning':'This inspects a frozen W90Flow package only; no NMEE response computed.'}

def main():
    import argparse
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('package');ap.add_argument('--structure');ap.add_argument('--output',type=Path)
    opts=ap.parse_args()
    result=inspect_package(opts.package,structure=opts.structure)
    body=json.dumps(result,indent=2,ensure_ascii=False)
    if opts.output:
        opts.output.parent.mkdir(parents=True,exist_ok=True);opts.output.write_text(body+'\n',encoding='utf-8')
    print(body)

if __name__=='__main__':main()
