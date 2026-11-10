"""Reanalyse all existing layer-channel controls; no new electronic calculation."""
from pathlib import Path
import json,csv,hashlib,argparse
import numpy as np
ROOT=Path(__file__).resolve().parent
parser=argparse.ArgumentParser();parser.add_argument('--materials-root',type=Path,default=ROOT.parent)
args=parser.parse_args()
specs={'MoSSe':('mosse-screened-channels-20261109/result.json','fm010','fp010',112),'InSe':('inse-nmee-20261109/channels.json','Am','Ap',160)}
rows=[];sources=[];checks={};paired=[]
for mat,(path,minus,plus,nbmain) in specs.items():
 src=args.materials_root/path;d=json.loads(src.read_text());assert d['complete']
 sources.append({'material':mat,'path':path,'sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'original_job':d['job'],'original_node':d['node'],'input_rows':len(d['rows'])})
 data=d['rows'];keys=sorted({(r['nb'],r['frozen_empty'],r['channel'],r['alpha']) for r in data})
 maxsum=0.
 for nb,frozen,ch,alpha in keys:
  rs={r['case']:r for r in data if (r['nb'],r['frozen_empty'],r['channel'],r['alpha'])==(nb,frozen,ch,alpha)}
  assert len(rs)==3 and all(k in rs for k in ('zero',minus,plus))
  cm=np.array(rs[minus]['C_total_lower_upper']);cp=np.array(rs[plus]['C_total_lower_upper']);c0=np.array(rs['zero']['C_total_lower_upper']);g=(cp-cm)/.2
  amp=(abs(cp)+abs(cm))/2
  assert np.all(amp>0)
  for c in (cm,cp,c0,g):maxsum=max(maxsum,abs(c[0]-c[1]-c[2]))
  r={'material':mat,'nb':nb,'frozen_empty':frozen,'channel':ch,'alpha':alpha,'representative':nb==nbmain and frozen==8,
     'minimum_gap_eV':min(v['gap_eV'] for v in rs.values())}
  for i,l in enumerate(('net','lower','upper')):
   r.update({f'C0_{l}':float(c0[i]),f'Cminus_{l}':float(cm[i]),f'Cplus_{l}':float(cp[i]),f'G_{l}':float(g[i]),
      f'Amean_{l}':float(amp[i]),f'R_{l}':float(abs(g[i])/amp[i]),
      f'amplitude_change_percent_{l}':float(100*(abs(cp[i])-abs(cm[i]))/amp[i]),
      f'endpoint_sign_change_{l}':bool(cp[i]*cm[i]<=0)})
  r['upper_share_abs_G']=float(abs(g[2])/(abs(g[1])+abs(g[2])))
  r['relative_upper_lower_ratio']=r['R_upper']/r['R_lower']
  rows.append(r)
 checks[mat]={'max_layer_sum_error':maxsum,'min_gap_eV':min(r['gap_eV'] for r in data),
   'max_half_probe_C_change':max(r.get('half_probe_max_change',0) for r in data),
   'max_transport_deviation':max(r['representation_max_singular_deviation'] for r in d['checks']),
   'row_count':len(keys),'all_windows':sorted({(r['nb'],r['frozen_empty']) for r in data})}
 assert maxsum<1e-6
for r in rows:
 full=next(x for x in rows if (x['material'],x['nb'],x['frozen_empty'],x['channel'])==(r['material'],r['nb'],r['frozen_empty'],'full'))
 for l in ('lower','upper'):
  r[f'abs_G_ratio_to_full_{l}']=abs(r[f'G_{l}']/full[f'G_{l}'])
  r[f'R_ratio_to_full_{l}']=r[f'R_{l}']/full[f'R_{l}']
  r[f'R_fixed_full_amplitude_{l}']=abs(r[f'G_{l}'])/full[f'Amean_{l}']
  r[f'C0_sign_differs_from_full_{l}']=r[f'C0_{l}']*full[f'C0_{l}']<0
ranges=[]
for mat in specs:
 for ch,a in sorted({(r['channel'],r['alpha']) for r in rows if r['material']==mat}):
  subset=[r for r in rows if (r['material'],r['channel'],r['alpha'])==(mat,ch,a)]
  metrics=['C0_lower','C0_upper','G_lower','G_upper','R_lower','R_upper','amplitude_change_percent_lower','amplitude_change_percent_upper','relative_upper_lower_ratio','upper_share_abs_G','abs_G_ratio_to_full_upper','R_ratio_to_full_upper']
  ranges.append({'material':mat,'channel':ch,'alpha':a,'n_windows':len(subset),
    'ranges':{k:[min(r[k] for r in subset),max(r[k] for r in subset)] for k in metrics},
    'all_lower_zero_gate_sign_changes':all(r['C0_sign_differs_from_full_lower'] for r in subset),
    'any_endpoint_layer_sign_change':any(r[f'endpoint_sign_change_{l}'] for r in subset for l in ('lower','upper'))})
for filename,rr in [('all_controls.csv',rows),('representative.csv',[r for r in rows if r['representative']])]:
 with (ROOT/filename).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rr[0]),lineterminator='\n');w.writeheader();w.writerows(rr)
out={'new_DFT':False,'new_remote_jobs':False,'sources':sources,'checks':checks,'ranges':ranges,'rows':rows,
 'definitions':{'C':'1e-9 (nm V T)^-1, static Ex Bz spin response','G':'[C(+.1)-C(-.1)]/.2; 1e-9 (V^2 T)^-1','Amean':'(|Cplus|+|Cminus|)/2','R':'|G|/Amean, (V/nm)^-1','amplitude_change_percent':'100*(|Cplus|-|Cminus|)/Amean; not change from zero gate'}}
(ROOT/'summary.json').write_text(json.dumps(out,indent=2)+'\n')
print('All',len(rows),'control/window groups analysed; source SCFs unchanged.')
for r in rows:
 if r['representative'] and r['channel'] in ('full','VC','Vb_Ct') and r['alpha'] in (0,1):
  print(r['material'],r['channel'],'C0_lower',round(r['C0_lower'],3),'upper |G| / full',round(r['abs_G_ratio_to_full_upper'],3),'upper R / full',round(r['R_ratio_to_full_upper'],3),'lower/upper amplitude %',round(r['amplitude_change_percent_lower'],3),round(r['amplitude_change_percent_upper'],3))
