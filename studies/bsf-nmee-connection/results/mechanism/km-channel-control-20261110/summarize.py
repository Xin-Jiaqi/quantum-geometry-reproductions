from pathlib import Path
import json,csv,hashlib,sys
import numpy as np
ROOT=Path(__file__).resolve().parent
files=['survey.json','refined.json','restoration.json','weak_training.json','weak_validation.json','phase_validation.json','contrast_validation.json']
rows=[];sources=[]
for f in files:
 d=json.loads((ROOT/f).read_text());assert d['complete'];sources.append({'file':f,'sha256':hashlib.sha256((ROOT/f).read_bytes()).hexdigest(),'grid':d['parameters']['grid'],'b':d['parameters']['b'],'du_meV':d['parameters']['du']})
 for r in d['rows']:
  row={k:r[k] for k in ['t_meV','U_meV','channel','retained_alpha','common_gap_eV','insulating','nk']};row['source']=f
  for i,l in enumerate(['lower','upper','net']):row['C_'+l]=r['C0'][i];row['G_'+l]=r['G'][i]
  rows.append(row)
with (ROOT/'all_results.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
pred=json.loads((ROOT/'weak_prediction.json').read_text());predmap={(r['t_meV'],r['alpha']):np.array(r['G']) for r in pred['predictions']}
errs=[]
for r in rows:
 if r['source']!='weak_validation.json' or r['t_meV']==0:continue
 expect=predmap[r['t_meV'],r['retained_alpha']];actual=np.array([r['G_'+l] for l in ['lower','upper','net']])
 errs.append({'t_meV':r['t_meV'],'alpha':r['retained_alpha'],'prediction':expect.tolist(),'actual':actual.tolist(),'relative_error':(abs(expect-actual)/abs(actual)).tolist(),'signs_correct':bool(np.all(expect*actual>0))})
ph=json.loads((ROOT/'phase_prediction.json').read_text());r=next(r for r in rows if r['source']=='phase_validation.json' and r['channel']=='Vb_Ct');pv=np.array(ph['G_predicted']);av=np.array([r['G_'+l] for l in ['lower','upper','net']]);phaseerr=(abs(pv-av)/abs(av)).tolist()
contrast=[]
for t in [10,20]:
 rr=[r for r in rows if r['source']=='contrast_validation.json' and r['t_meV']==t];full=next(r for r in rr if r['channel']=='full');off=next(r for r in rr if r['channel']=='Vb_Ct')
 contrast.append({'t_meV':t,'U_meV':-4,'G_upper_full':full['G_upper'],'G_upper_channel_removed':off['G_upper'],'full_to_removed_abs_G_ratio':abs(full['G_upper']/off['G_upper'])})
summary={'complete':True,'sources':sources,'rows':rows,'weak_validation':errs,'phase_relative_error':phaseerr,'phase_signs_correct':bool(np.all(pv*av>0)),'strong_contrast':contrast,'minimum_gap_eV':min(r['common_gap_eV'] for r in rows),'all_insulating':all(r['insulating'] for r in rows),'layer_sum_C_max':max(abs(r['C_lower']+r['C_upper']-r['C_net']) for r in rows),'layer_sum_G_max':max(abs(r['G_lower']+r['G_upper']-r['G_net']) for r in rows),'max_weak_layer_relative_error':max(max(e['relative_error'][:2]) for e in errs),'scope':'No DFT, no added states. Relative sign and channel scaling are diagnostic Hamiltonian interventions, not experimentally validated knobs.'}
(ROOT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print('All rows',len(rows),'gap',summary['minimum_gap_eV'],'max weak layer error',summary['max_weak_layer_relative_error']);print('contrast',contrast);print('phase sign',summary['phase_signs_correct'])
