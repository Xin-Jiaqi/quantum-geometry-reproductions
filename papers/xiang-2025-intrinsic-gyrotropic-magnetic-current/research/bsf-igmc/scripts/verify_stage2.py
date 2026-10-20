"""Offline, clean-module acceptance: tests plus locked two-branch benchmark."""
import subprocess,sys,json,hashlib,os
from pathlib import Path
base=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(base/'src'))
from bsf_igmc_micro.model import Parameters,geometric_response,layer_polarization,connection_error
from bsf_igmc_micro.symmetry import screen

def verify():
  lock=json.loads((base/'results/PROVENANCE.json').read_text())
  for path,h in lock['sha256'].items():
    data=(base/path).read_bytes()
    actual=hashlib.sha256(data).hexdigest()
    if actual!=h:raise ValueError(f'SOURCE_DRIFT {path}')
  for short,branch,mu,sign in [('af','AF',0.,-1),('fm','FM',-.03,1)]:
    data=json.loads((base/f'results/stage2_{short}.json').read_text())['branches'][branch]
    model=Parameters(**json.loads((base/f'results/stage2_{short}.json').read_text())['parameters'])
    recalc=geometric_response(.19,branch,model,N=21,mu=mu)
    ref=data['plus']['21']
    if abs(recalc['Dxx_over_omega']-ref['Dxx_over_omega'])>5e-11:raise ValueError('NUMERIC_DRIFT '+branch)
    if recalc['indirect_gap']<=0 or recalc['occ_range']!=[4,4]:raise ValueError('NO_GAP '+branch)
    other=data['minus']['31']['Dxx_over_omega']; ref31=data['plus']['31']['Dxx_over_omega']
    if abs(ref31-sign*other)>1e-10:raise ValueError('PARTNER_PARITY '+branch)
    if abs(data['layer_charge_+0.19']+data['layer_charge_-0.19'])>1e-10:raise ValueError('DIPOLE_PARITY '+branch)
  env=os.environ.copy();env['PYTHONPATH']=str(base/'src')+os.pathsep+env.get('PYTHONPATH','')
  env['OPENBLAS_NUM_THREADS']='1'
  tests=subprocess.run([sys.executable,'-m','unittest','discover','-s','tests','-v'],cwd=base,env=env,check=True)
  return {'status':'PASS_MODEL_LEVEL_ONLY','physics_status':'BSF_FE_BLOCKED_REAL_MATERIAL_BLOCKED','tests':'all dynamically discovered unittest tests passed'}

if __name__=='__main__':print(json.dumps(verify(),indent=2))
