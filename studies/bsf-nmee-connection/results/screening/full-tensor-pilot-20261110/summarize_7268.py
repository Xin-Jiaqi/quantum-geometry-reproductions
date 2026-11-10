"""Summarize returned small result JSONs; no DFT/cache processing or dispatch."""
from pathlib import Path
import json
import numpy as np
ROOT=Path(__file__).resolve().parent

def zero_check():
    contract=json.loads((ROOT/'zero_acceptance_contract.json').read_text())
    rows=[]
    for case in ('A0','lower0'):
        p=ROOT/'pilot7268'/case/'result.json'
        if not p.exists():return {'passed':False,'status':'waiting_for_results','missing':case}
        d=json.loads(p.read_text())
        if not d.get('complete'):return {'passed':False,'status':'incomplete','case':case}
        scopes=[0,1,2] if case=='A0' else [0]
        checks=[];signs=[]
        for nb,res in d['responses'].items():
            r=res['standard'];c=np.array(r['C'])[scopes];half=np.array(res['half_b']['C'])[scopes]
            scale=max(float(np.max(abs(c[:,0,0]))),1e-8)
            forbidden=float(np.max(abs(c[:,[0,1,2,2],[1,0,0,1]])))
            iso=float(np.max(abs(c[:,0,0]-c[:,1,1])))
            probe=float(np.max(abs(c-half)))
            checks.append({'nb':nb,'Cxx_total_lower_upper' if case=='A0' else 'Cxx_isolated':c[:,0,0].tolist(),'scale':scale,'forbidden_max_abs':forbidden,'Cxx_Cyy_max_abs':iso,'probe_max_abs':probe,'gap_eV':r['gap_eV'],'passed':bool(r['gap_eV']>contract['min_gap_eV'] and forbidden<.02*scale and iso<.05*scale and probe<.01*scale)})
            signs.append(np.sign(c[:,0,0]))
        same_signs=all(np.array_equal(signs[0],x) for x in signs[1:])
        force_pass=d['max_force_eV_A']<=contract['max_force_eV_A']
        rows.append({'case':case,'checks':checks,'same_signs_across_band_windows':same_signs,'max_force_eV_A':d['max_force_eV_A'],'force_pass':force_pass,'passed':all(x['passed'] for x in checks) and same_signs and force_pass,'scf_seconds':d['scf_seconds'],'total_seconds':d['elapsed_seconds']})
    return {'passed':all(x['passed'] for x in rows),'status':'reviewed','rows':rows,'upper_reference':'same as lower by verified pure in-plane translation','scope':'fixed raw database-derived geometry, full shell method; not final mesh convergence'}

if __name__=='__main__':
    result=zero_check();(ROOT/'zero7268_acceptance.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
