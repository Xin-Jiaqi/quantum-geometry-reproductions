"""Explicit planning scenarios, not measured runtime claims for uncomputed materials."""
from pathlib import Path
import csv,json,math
import numpy as np
from ase.io import read
ROOT=Path(__file__).resolve().parent
PLAN=ROOT.parent/'bidb-pilot-plan-20261110'
records=json.loads((ROOT/'historical_scf_timings.json').read_text())['records']
def median(material,kind):return float(np.median([r['scf_seconds'] for r in records if r['material']==material and r['kind']==kind]))
ref=read(PLAN/'structures/bilayer-7269.vasp');refV=ref.get_volume();spacing=max(np.linalg.norm(2*np.pi*np.linalg.inv(ref.cell.array).T,axis=1)[:2])/6
# Electron counts are planning proxies. Exact licensed PAW sets must be fixed at dispatch.
valence={'Mo':12,'W':12,'In':13,'Ga':13,'Ge':4,'Sn':14,'S':6,'Se':6}
rows=[]
for candidate in csv.DictReader((PLAN/'candidates.csv').open()):
 uid=int(candidate['id']);at=read(PLAN/f'structures/bilayer-{uid}.vasp')
 rec=2*np.pi*np.linalg.inv(at.cell.array).T
 mesh=np.maximum(4,np.ceil(np.linalg.norm(rec,axis=1)[:2]/spacing-1e-8).astype(int));nk=int(np.prod(mesh))
 ne=sum(valence[s] for s in at.get_chemical_symbols())
 volratio=at.get_volume()/refV
 factor=nk/36*volratio*(ne/48)**2
 material=candidate['parent_formula'].replace('In2Se2','InSe').replace('Ga2Se2','GaSe').replace('Ge2S2','GeS').replace('Sn2S2','SnS')
 anchor='InSe' if material in ('InSe','GaSe') else 'MoSSe'
 if anchor=='InSe':
  inse=read(PLAN/'structures/bilayer-864.vasp')
  factor=nk/36*at.get_volume()/inse.get_volume()*(ne/76)**2
  tb,tm=median('InSe','bilayer')*factor,median('InSe','mono')*factor;reason='InSe family historical timings scaled by Nk*cell volume*electron-proxy squared'
 else:
  tb,tm=median('MoSSe','bilayer')*factor,median('MoSSe','mono')*factor;reason='MoSSe timing scaled by Nk*cell volume*electron-proxy squared'
 # Cold zero-field iteration overhead is unmeasured for most new candidates.
 # Nominal 1.5x only for the 3 zero-field states; two gates retain anchor medians.
 post_bilayer=35*nk/36*(ne/48)**3
 post_mono=35*nk/36*((ne/2)/48)**3
 nominal=(3.5*(tb+2*tm)+3*post_bilayer+6*post_mono)/3600
 lowfactor=.5;highfactor=2.
 if material in ('GeS','SnS'):
  lower_ne=sum((4 if s in ('Ge','Sn') else valence[s]) for s in at.get_chemical_symbols())
  upper_ne=sum((14 if s in ('Ge','Sn') else valence[s]) for s in at.get_chemical_symbols())
  lowfactor*= (lower_ne/ne)**2;highfactor*= (upper_ne/ne)**2
 rows.append({'id':uid,'material':material,'mesh_planning_only':f'{mesh[0]}x{mesh[1]}','states_no_reuse':9,'nominal_serial_job_hours':nominal,'planning_low_hours':nominal*lowfactor,'planning_high_hours':nominal*highfactor,'nominal_allocated_CPU_hours_4rank':nominal*4,'status':'scenario_not_precision_prediction','basis':reason})
with (ROOT/'runtime_estimates.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
summary={'scope':'16 bilayers plus2geometry-matched isolated references each, 3F states each =144 SCF states before any reuse; fixed ions, exploratory reciprocal spacing, no relaxation/convergence campaign',
 'time_definition':'sum of job elapsed pure compute hours at4CPU; excludes queue/transfer; CPU-hour equals4times elapsed allocated walltime, not measured CPU utilization',
 'model':'anchored in historical MoSSe/InSe SCF times and actual7269 tensor processing; unseen systems scale Nk*volume*electron-count-proxy^2; cold zero-SCF factor1.5; each postprocessing state scales Nk*its own electron-count^3 (3 bilayer and6 monolayer states); deliberately broad0.5–2 envelope, plus Ge/Sn PAW ambiguity; not confidence intervals',
 'limits':['No benchmark for6new chemical families','Final actual PAW sets and Nb not frozen','GeS/SnS symmetry and grid acceptance unresolved','No guaranteed convergence at planned mesh','No savings from cached SCFs/reference dedup assumed','No ionic relaxation, NEB, dense grids or mechanism intervention budget','Two concurrent4CPU jobs give ideal elapsed lower bound sum/2, not schedule guarantee'],
 'totals':{k:sum(r[k] for r in rows) for k in ['states_no_reuse','nominal_serial_job_hours','planning_low_hours','planning_high_hours','nominal_allocated_CPU_hours_4rank']},
 'by_material':{m:{'bilayers':len([r for r in rows if r['material']==m]),'per_bilayer_nominal_hours':float(np.mean([r['nominal_serial_job_hours'] for r in rows if r['material']==m])),'per_bilayer_low_hours':min(r['planning_low_hours'] for r in rows if r['material']==m),'per_bilayer_high_hours':max(r['planning_high_hours'] for r in rows if r['material']==m)} for m in dict.fromkeys(r['material'] for r in rows)}}
(ROOT/'runtime_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(summary['by_material'],ensure_ascii=False,indent=2));print(summary['totals'])
