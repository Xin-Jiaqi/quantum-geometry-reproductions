#!/usr/bin/env python3
"""Strict source-versus-model comparison of *printed* SM Tables I/II, pp8-10.

Runs independently of the author's full group theory toolkit installation.
Return code 0 if 243 expected matches and the *one known printed-row
independence omission* is reproduced exactly. Any new discrepancy fails.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from igmc_lab.paper_table_source import (ALL_SOURCE_ROWS,C_SOURCE_ROWS,D_SOURCE_ROWS,
    compare_source_table,source_in_toolkit_setting)
from igmc_lab.toolkit_bridge import read_toolkit_snapshot,tensor_basis

KNOWN = ('D',"-42'm'")

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--no-write',action='store_true')
    a=p.parse_args(argv)
    groups=read_toolkit_snapshot()
    records,summary=compare_source_table(groups)
    mismatches=summary['mismatches']
    if len(mismatches)!=1 or (mismatches[0]['channel'],mismatches[0]['symbol'])!=KNOWN:
        raise AssertionError(f'Unexpected mismatches: {mismatches}')
    m=mismatches[0]
    if m['sm_dimension']!=2 or m['independently_calculated_dimension']!=1:
        raise AssertionError('Known exception changed; do not silently accept.')
    # Evaluate whether the source row has right ALLOWED component support,
    # although it declares two independently-variable xy/yx parameters.
    group=next(x for x in groups if x.symbol==KNOWN[1]);base=tensor_basis(group,'D')
    sheet=next(row for row in D_SOURCE_ROWS if row.id=='D05')
    A=source_in_toolkit_setting(sheet)
    # A 1-dim permitted subspace must be fully contained in the printed 2-dim one.
    qa,_=np.linalg.qr(A)
    contained_error=float(np.max(abs(base-(qa@qa.T@base))))
    if contained_error>1e-10:raise AssertionError('Source pattern not even a superset')
    # Self-contained analytic independent check using the explicit generators:
    # for unprimed S4, det=-1 and R=[[0,1,0],[-1,0,0],[0,0,-1]]
    # -> planar T=[[a,b],[c,d]] maps to [[-d,c],[b,-a]], b=c.
    # for primed C2x, eta_D=-1 and R=diag(1,-1,-1)
    # -> planar T maps to [[-a,b],[c,-d]], forcing a=d=0.
    # Together: Txy=Tyx, all others zero.
    symmetry=np.zeros((9,1));symmetry[1,0]=1;symmetry[3,0]=1
    q,_=np.linalg.qr(base)
    analytic_residual=float(np.max(abs((q@q.T)@symmetry-symmetry)))
    if analytic_residual>1e-10:
        raise AssertionError('Analytic S4/C2x constraints disagree with MPG data')
    summary['exception']={
        'channel':KNOWN[0],'group':KNOWN[1],'magnetic_group_number':group.number,
        'sm_table':'II','sm_page':9,'sm_row':sheet.id,
        'printed_pattern':'sigma_xy and sigma_yx listed as independent',
        'calculated_pattern':'sigma_xy = sigma_yx; all other components vanish',
        'permitted_tensor_shape_preserved':True,
        'printed_dimension':2,'symmetry_dimension':1,
        'source_superset_containment_error':contained_error,
        'independent_generator_analytic_subspace_residual':analytic_residual,
        'classification_affected':False,
        'conclusion':'Missing equality relation in printed Table II row; do not overstate as nonzero/forbidden mismatch.',
    }
    summary['status']='PASS_WITH_ONE_SOURCE_TABLE_INDEPENDENCE_EXCEPTION'
    summary['source_path']='supple-2(1).pdf, Table I pp8-9 / Table II pp9-10; manually transcribed and visually checked'
    summary['source_transcription_sha256']=hashlib.sha256((ROOT/'src/igmc_lab/paper_table_source.py').read_bytes()).hexdigest()
    if not a.no_write:
        root=ROOT/'results/raw';root.mkdir(exist_ok=True,parents=True)
        with (root/'round5_source_table_244.csv').open('w',newline='',encoding='utf8') as f:
            writer=csv.DictWriter(f,fieldnames=records[0]);writer.writeheader();writer.writerows(records)
        (root/'round5_source_table_audit.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
    print(json.dumps(summary,indent=2,ensure_ascii=False))
    return 0
if __name__=='__main__':raise SystemExit(main())
