"""Reevaluate the archived gate checks and the frozen prediction from accepted records."""
import csv
import json
from pathlib import Path
import sys
import numpy as np

here = Path(__file__).resolve().parent
sys.path.insert(0,str(here.parent/'overnight7268-20261110'))
from night_checks import response_acceptance

def save(name, obj):
    (here/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')

records = {}
table = []
for path in sorted((here/'states').glob('*/result.json')):
    row = json.loads(path.read_text())
    limit = .03 if row['label'] == 'A' and row['Fz_V_nm'] == 0 else .1
    accepted = response_acceptance(row, limit)
    assert accepted['passed'] and row['complete'] and row['acceptance']['passed']
    records[row['case']] = row
    for nb, response in row['responses'].items():
        c = np.array(response['standard']['C'])
        scopes = ['bilayer_total','bilayer_lower','bilayer_upper'] if row['label']=='A' else ['isolated_'+row['label']]
        for i,scope in enumerate(scopes):
            table.append(dict(mesh=row['mesh'],scope=scope,field_V_nm=row['Fz_V_nm'],bands=int(nb),
                              Cxx=c[i,0,0],Cyy=c[i,1,1],Cxy=c[i,0,1],Cyx=c[i,1,0],
                              Czx=c[i,2,0],Czy=c[i,2,1],gap_eV=response['standard']['gap_eV']))
with (here/'field_responses.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=table[0],lineterminator='\n');w.writeheader();w.writerows(table)

def tensor(mesh,label,suffix,nb,probe='standard'):
    return np.asarray(records[f'k{mesh}_{label}{suffix}']['responses'][str(nb)][probe]['C'])

metrics = []
for mesh in (6,9,12):
    for nb in (96,112):
        cm,c,cp = (tensor(mesh,'A',s,nb) for s in ('m','0','p'))
        g=(cp-cm)/.2
        gh=(tensor(mesh,'A','p',nb,'half_b')-tensor(mesh,'A','m',nb,'half_b'))/.2
        gl,gu=g[1:,0,0]
        metrics.append(dict(mesh=mesh,bands=nb,C0=c[:,0,0].tolist(),G=g[:,0,0].tolist(),
                            static_upper_share=float(abs(c[2,0,0])/sum(abs(c[1:,0,0]))),
                            gate_upper_share=float(abs(gu)/(abs(gu)+abs(gl))),
                            upper_minus_lower_slope=float(gu-gl),
                            probe_G_change=float(np.max(abs(gh[1:,0,0]-g[1:,0,0])))))
diag={}
for old,new in [(6,9),(9,12)]:
    row={}
    for q in ('C0','G'):
        def get(mesh,nb):
            if q=='C0':return tensor(mesh,'A','0',nb)[1:]
            return (tensor(mesh,'A','p',nb)-tensor(mesh,'A','m',nb))[1:]/.2
        x,y=get(old,112),get(new,112)
        delta=float(np.max(abs(y-x)));scale=float(np.max(abs(y[:,0,0])))
        band=float(max(np.max(abs(get(m,96)[:,0,0]-get(m,112)[:,0,0])) for m in (old,new)))
        margin=float(abs(abs(y[1,0,0])-abs(y[0,0,0])))
        if q=='G': probe=max(m['probe_G_change'] for m in metrics if m['mesh'] in (old,new))
        else:probe=max(float(np.max(abs(tensor(m,'A','0',n)[1:,0,0]-tensor(m,'A','0',n,'half_b')[1:,0,0]))) for m in (old,new) for n in (96,112))
        row[q]=dict(mesh_change=delta,relative=delta/scale,pilot_amplitude_pass=delta/scale<=(.05 if q=='C0' else .10),
                    band_change=band,probe_change=probe,rank_margin=margin,
                    rank_resolved=bool(margin>2*max(delta,band,probe) and all(abs(get(m,n)[1,0,0])>abs(get(m,n)[0,0,0]) for m in (old,new) for n in (96,112))))
    diag[f'{old}_to_{new}']=row
rows=[r for r in records.values() if r['label']=='A']
gates=[a for r in rows for a in r['acceptance']['rows']]
save('prediction_review.json',dict(source_job='565092.cmspbs1',accepted_state_count=len(records),
      raw_tensor_csv_units='1e-9 (nm V T)^-1; G in 1e-9 (V^2 T)^-1',
      mesh_metrics=metrics,bilayer_convergence=diag,
      min_bilayer_probe_gap_eV=min(a['gap_eV'] for a in gates),
      max_bilayer_force_eV_A=max(r['max_force_eV_A'] for r in rows),
      max_forbidden_per_scale=max(a['forbidden']/a['scale'] for a in gates),
      max_isotropy_per_scale=max(a['isotropy']/a['scale'] for a in gates),
      max_probe_change_per_scale=max(a['probe_difference']/a['scale'] for a in gates),
      max_layer_closure=max(a['layer_closure'] for a in gates),
      frozen_prediction='parallel-idea-loop-20261110/physics/prediction_freeze.json',
      verdict={'zero_directions':'supported','static_upper':'supported','gate_lower':'contradicted',
               'different_static_and_sensitive_layer':'contradicted','no_reversal_at_sampled_fields':'supported',
               'unsampled_fields':'no continuous field scan','positive_gate_amplitude_direction':'not predicted; both weaken'},
      missing_from_night=['k12_upperm','k12_upperp']))

ss_dir = here/'ss_comparison'
if (ss_dir/'fp010/result.json').exists():
    zero = json.loads((here.parent/'full-tensor-pilot-20261110/stencil_diagnosis.json').read_text())
    assert zero['passed']
    controls = {s: json.loads((ss_dir/s/'result.json').read_text()) for s in ('fm010','fp010')}
    for record in controls.values():
        assert response_acceptance(record,.1)['passed']
    compare = []
    for nb in (96,112):
        ss0 = next(row for row in zero['rows'] if row['nb']==nb)
        for interface in ('S/S','S/Se'):
            if interface=='S/S':
                a,b = [np.array(controls[s]['responses'][str(nb)]['standard']['C']) for s in ('fm010','fp010')]
                c = np.array(ss0['standard']['C'])
            else:
                a,c,b=[tensor(6,'A',s,nb) for s in ('m','0','p')]
            g=((b-a)/.2)[:,0,0]
            compare.append(dict(interface=interface,mesh=6,bands=nb,Cminus=a[:,0,0].tolist(),Czero=c[:,0,0].tolist(),
                                Cplus=b[:,0,0].tolist(),G=g.tolist(),upper_absolute_share=float(abs(g[2])/sum(abs(g[1:])))))
    save('matched_kernel_interfaces.json',compare)
