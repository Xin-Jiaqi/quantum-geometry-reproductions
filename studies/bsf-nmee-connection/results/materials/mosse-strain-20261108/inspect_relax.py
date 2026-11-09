from pathlib import Path
import os,socket,re,json,datetime
import numpy as np
if not os.getenv('PBS_JOBID') or socket.gethostname().startswith(('super','cms-ssh')):raise SystemExit('PBS compute only')
p=Path('ep02/relax');text=(p/'OUTCAR').read_text(errors='replace');part=text.rsplit('TOTAL-FORCE (eV/Angst)',1)[1];rows=[]
for line in part.splitlines():
 values=line.split()
 if len(values)==6:
  try:rows.append(list(map(float,values)))
  except ValueError:pass
 if len(rows)==6:break
assert len(rows)==6
arr=np.array(rows);z=np.sort(arr[:,2]);ionic=re.findall(r'^\s*(\d+) F=.*$',(p/'OSZICAR').read_text(),re.M)
out={'job':os.environ['PBS_JOBID'],'node':socket.gethostname(),'time':datetime.datetime.now().astimezone().isoformat(),'completed_ionic_steps':int(ionic[-1]),'positions_forces':rows,'max_z_force_eVA':float(abs(arr[:,5]).max()),'max_xy_force_eVA':float(abs(arr[:,3:5]).max()),'net_z_force_eVA':float(arr[:,5].sum()),'interlayer_gap_A':float(z[3]-z[2]),'ionic_converged':'reached required accuracy' in text}
Path('relax_inspection.json').write_text(json.dumps(out,indent=2));print(json.dumps(out),flush=True)
