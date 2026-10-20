#!/usr/bin/env python3
"""Cross-check 122 standard MPGs of Xiang SM Tables I–III with Xin-Jiaqi toolkit generators."""
from __future__ import annotations
import argparse,csv,hashlib,json,sys
from collections import Counter
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from igmc_lab.toolkit_bridge import (read_toolkit_snapshot,closure_order,tensor_basis,
                                     classify_basis,group_invariance_max,toolkit_live_dimensions,
                                     GENERATOR_FILE,canonical_component_basis,character_dimension)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--require-live-toolkit',action='store_true',help='also verify against an installed matching external group toolkit; fail if missing')
    args=ap.parse_args()
    gs=read_toolkit_snapshot();rows=[];tensor_shapes=[]
    for g in gs:
        order=closure_order(g)
        if order!=g.order:raise AssertionError('magnetic group fails closure',g.symbol,order,g.order)
        bC=tensor_basis(g,'C');bD=tensor_basis(g,'D')
        c_char=character_dimension(g,'C');d_char=character_dimension(g,'D')
        if c_char!=bC.shape[1] or d_char!=bD.shape[1]:
            raise AssertionError('SVD / character mismatch',g.symbol,c_char,d_char)
        row={'number':g.number,'symbol':g.symbol,'family':g.family,'order':order,
             'C_dim':bC.shape[1],'C_char_dim':c_char,'C_class':classify_basis(bC),
             'D_dim':bD.shape[1],'D_char_dim':d_char,'D_class':classify_basis(bD),
             'Eq14_as_C_dim':tensor_basis(g,'C',protocol='paper_printed_eq14').shape[1]}
        for kind,b in (('C',bC),('D',bD)):
            residual=group_invariance_max(g,b,kind)
            if residual>1e-6:raise AssertionError('tensor fails Neumann test',row,residual)
        rows.append(row)
        shape={'number':g.number,'symbol':g.symbol}
        for kind,b in (('C',bC),('D',bD)):
            canonical=canonical_component_basis(b)
            if group_invariance_max(g,canonical,kind)>1e-6:
                raise AssertionError('canonical basis broken',g.symbol,kind)
            shape[kind+'_basis']=[[[float(round(x,7)) for x in rr] for rr in canonical[:,i].reshape(3,3)]
                                   for i in range(canonical.shape[1])]
        tensor_shapes.append(shape)
    counts={k:dict(Counter(row[f'{k}_class'] for row in rows)) for k in ('C','D')}
    assert counts['C']=={'both':16,'transverse_only':18,'longitudinal_only':24,'forbidden':64},counts['C']
    assert counts['D']=={'both':16,'transverse_only':19,'longitudinal_only':23,'forbidden':64},counts['D']
    live='not_executed_external_toolkit_unavailable'
    if args.require_live_toolkit:
        differences=toolkit_live_dimensions(gs)
        if differences:raise AssertionError(differences)
        live='PASS_244_external_toolkit_API_dimensions'
    raw=ROOT/'results'/'raw'
    raw.mkdir(parents=True,exist_ok=True)
    columns=['number','symbol','family','order','C_dim','C_char_dim','C_class','D_dim','D_char_dim','D_class','Eq14_as_C_dim']
    with (raw/'round4_mpg_122.csv').open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=columns);writer.writeheader();writer.writerows(rows)
    (raw/'round4_mpg_tensor_shapes.json').write_text(json.dumps(tensor_shapes,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    differences=[{'number':r['number'],'symbol':r['symbol'],'C_dim':r['C_dim'],'Eq14_C_dim':r['Eq14_as_C_dim']}
                 for r in rows if r['C_dim']!=r['Eq14_as_C_dim']]
    report={'provenance':{'github_repo':'Xin-Jiaqi/group-theory-operations-toolkit',
            'git_tree':'551a0daf8d4ba5a5f24b2811694bfe438616de48',
            'magnetic_point_groups_blob_sha':'406309c2686cbdaf7634dbea10ad8b059ee40e4d',
            'group_operations_blob_sha':'1c43c66646e214c8f5939c06a17a8e69fb4c5122',
            'minimal_snapshot_sha256':hashlib.sha256(GENERATOR_FILE.read_bytes()).hexdigest()},
            'protocol':'real AC in-phase conduction C / quadrature displacement D; not full complex finite-omega tensor',
            'method':'122 complete magnetic group closures; 244 independent 9D Neumann-SVD subspaces checked against 244 complete-group character averages',
            'status':'PASS_METHOD_SYMMETRY_SUBSET','character_crosschecks':244,'class_counts':counts,'allowed_groups':{'C':58,'D':58},
            'Eq14_C_dimension_disagreements':len(differences),
            'Eq14_mismatches':differences,
            'matched_sm_reference':{'4mm1prime':{'C_dim':1,'D_dim':0},'4mprimemprime':{'C_dim':1,'D_dim':2}},
            'external_live_toolkit_comparison':live,
            'coordinate_warning':'SM monoclinic 2/m uses unique axis y; standard toolkit snapshots use z; shape comparison requires y<->z transform'}
    (raw/'round4_mpg_audit.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    reference={'toolkit_snapshot_sha256':report['provenance']['minimal_snapshot_sha256'],
               'csv_sha256':hashlib.sha256((raw/'round4_mpg_122.csv').read_bytes()).hexdigest(),
               'tensor_shapes_sha256':hashlib.sha256((raw/'round4_mpg_tensor_shapes.json').read_bytes()).hexdigest(),
               'class_counts':counts,'Eq14_C_dim_mismatch':len(differences),
               'method':'closure+Neumann+character average, 122 groups 244 tensors, NOT author plotting values'}
    frozen=ROOT/'results'/'frozen'/'round4_mpg_reference_v2.json'
    if frozen.exists():
        if json.loads(frozen.read_text(encoding='utf-8'))!=reference:raise AssertionError('Round 4 reference drift')
    else:frozen.write_text(json.dumps(reference,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='Eq14_mismatches'},indent=2,ensure_ascii=False))
if __name__=='__main__':main()
