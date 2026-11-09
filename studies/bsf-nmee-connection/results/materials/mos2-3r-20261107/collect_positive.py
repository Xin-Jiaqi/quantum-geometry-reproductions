"""Collect completed +/-E densities at positive B using nonmagnetic time reversal."""
import ast,json,os,socket,hashlib
from pathlib import Path
import numpy as np
if not os.getenv('PBS_JOBID') or not os.getenv('PBS_NODEFILE') or socket.gethostname().startswith(('super','cms-ssh')):raise SystemExit('compute node required')
root=Path.cwd();label=os.environ['STRUCTURE'];parent=os.environ['SOURCE_JOB'];gate=0.;electric=.004;zeeman=.016;mesh=6
base=root/f'{label}_F{gate:+.3f}_k{mesh}_e{electric}_b{zeeman}'
fn=next(n for n in ast.parse((root/'response.py').read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='spin')
ns={'np':np,'ne':48};exec(compile(ast.Module(body=[fn],type_ignores=[]),'spin','exec'),ns)
rows={};gaps={}
for name in ['E+1B+1','E-1B+1']:
 d=base/name;finished=json.loads((d/'finished.json').read_text())
 if finished['gap_eV']<=0:raise ValueError('gap closed')
 rows[name]=ns['spin'](d);gaps[name]=finished['gap_eV']
lines=(root/f'POSCAR.{label}').read_text().splitlines();cell=np.array([list(map(float,s.split())) for s in lines[2:5]])*float(lines[1]);area=float(np.linalg.norm(np.cross(cell[0],cell[1])))
mu=2.00231930436256*5.7883818060e-5/2;c=(np.array(rows['E+1B+1'])[1:]-np.array(rows['E-1B+1'])[1:])*(1e10*mu/area)/(2*electric*zeeman)
if np.max(abs(c[:,0]-c[:,1]-c[:,2]))>1e-7:raise ValueError('layer sum')
r={'structure':label,'gate_physical_Vnm':gate,'vasp_EFIELD':0.,'probe_VA':electric,'zeeman_eV':zeeman,'mesh':mesh,'area_A2':area,'C_xyz_total_bottom_top_1e-9':c.tolist(),'B_even_residual_xyz':None,'one_B_sign_C_xyz':[c.tolist()],'B_signs':[1],'spin_raw':rows,'gaps_eV':gaps,'pbs':parent,'analysis_pbs':os.environ['PBS_JOBID'],'node':socket.gethostname(),'structure_sha256':hashlib.sha256((root/f'POSCAR.{label}').read_bytes()).hexdigest(),'mono_halfspace_warning':False,'protocol_note':'Nonmagnetic time reversal; original job deliberately stopped after completed positive-B +/-E pair. No negative-B numerical residual claimed.'}
(base/'answer.json').write_text(json.dumps(r,indent=2));print('ANSWER',json.dumps(r))
