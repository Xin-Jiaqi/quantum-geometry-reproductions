"""Small-result integrity checks and physical comparisons; no new DFT."""
from pathlib import Path
import hashlib,json
import numpy as np
from ase.io import read
p=Path(__file__).resolve().parent
jobs=json.loads((p/'jobs.json').read_text());status=json.loads((p/'final_job_status.json').read_text())
baseline=json.loads((p/'baseline.json').read_text());data=[baseline];audits=[]
for case,job in jobs.items():
    d=p/case;r=json.loads((d/'result.json').read_text());e=json.loads((d/'export.json').read_text())
    assert status[job.split('.')[0]]['Exit_status']=='0'
    assert r['job']==e['job']==job and r['case']==case
    assert r['node'].startswith('csuper') and e['node'].startswith('csuper')
    for name,expected in r['input_sha256'].items():
        assert hashlib.sha256((d/'seed'/name).read_bytes()).hexdigest()==expected
    atoms=read(d/'seed/POSCAR');lo=atoms.positions[:,2]<atoms.cell[2,2]/2
    gap=float(atoms.positions[~lo,2].min()-atoms.positions[lo,2].max())
    assert abs(gap-r['interface_gap_A'])<1e-10
    assert abs((gap-r['delta_d_A'])-3.03829445)<1e-7
    inc=(d/'seed/INCAR').read_text();assert 'BEXT=0 0 0' in inc and 'EFIELD_PEAD' not in inc
    audit={'case':case,'step_change_max':0.,'gap_min_eV':e['gap_eV']}
    for nb,v in r['results'].items():
        x=np.array(v['standard']['C_total_lower_upper']);y=np.array(v['half_b']['C_total_lower_upper'])
        assert np.isfinite(x).all() and abs(x[0]-x[1]-x[2])<1e-9
        assert v['standard']['gap_eV']>0 and v['half_b']['gap_eV']>0
        audit['step_change_max']=max(audit['step_change_max'],float(abs(x-y).max()))
        audit['gap_min_eV']=min(audit['gap_min_eV'],v['standard']['gap_eV'],v['half_b']['gap_eV'])
    audits.append(audit);data.append(r)
rows=[]
for r in sorted(data,key=lambda q:(q['delta_d_A'],q['physical_Fz_Vnm'])):
    x=np.array([r['results'][nb]['standard']['C_total_lower_upper'] for nb in ('96','112')])
    rows.append({'case':r['case'],'delta_d_A':r['delta_d_A'],'Fz_Vnm':r['physical_Fz_Vnm'],
                 'C_96_total_lower_upper':x[0].tolist(),'C_112_total_lower_upper':x[1].tolist(),
                 'signs_agree_between_cutoffs':bool(np.all(np.sign(x[0])==np.sign(x[1])))})
comparisons={}
for nb in ('96','112'):
    lookup={(r['delta_d_A'],r['physical_Fz_Vnm']):np.array(r['results'][nb]['standard']['C_total_lower_upper']) for r in data}
    comparisons[nb]={
        'gate_change_at_d0':(lookup[0.,-.5]-lookup[0.,0.]).tolist(),
        'gate_change_at_d030':(lookup[.3,-.5]-lookup[.3,0.]).tolist(),
        'distance_change_at_zero_gate':(lookup[.3,0.]-lookup[0.,0.]).tolist(),
        'distance_change_at_minus_half_gate':(lookup[.3,-.5]-lookup[0.,-.5]).tolist()}
summary={'rows':rows,'audit':audits,'comparisons_total_lower_upper':comparisons,
         'units':'1e-9 (nm V T)^-1, Pauli spin density',
         'scope':'All actual points; no fitted critical distance/field. Positive gaps are checked at sampled points, not all interpolating structures.'}
(p/'summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
