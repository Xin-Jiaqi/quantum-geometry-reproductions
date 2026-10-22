"""v70 multi-frequency joint finite-lead-length / spectral-resolution audit."""
import argparse,json,os,subprocess,sys
from pathlib import Path
import numpy as np
from qgmap.finite_lead_lehmann_v69 import Spec,lehmann_response
from qgmap.lead_keldysh_eb_v68 import Spec as KSpec,response as keldysh_response

FREQUENCIES=(.20,.29,.40)

def production_value(omega, length=1000, eta=.005):
    """Run one large diagonalization in a child process to release LAPACK workspace."""
    root=Path(__file__).resolve().parents[1]
    code=("import json,sys; from qgmap.finite_lead_lehmann_v69 import Spec,lehmann_response; "
          "z=lehmann_response(Spec(omega=float(sys.argv[1]),length=int(sys.argv[2]),eta=float(sys.argv[3])))['sigma_device']; "
          "print(json.dumps({'real':float(z.real),'imag':float(z.imag)}))")
    env=os.environ.copy();env['OPENBLAS_NUM_THREADS']='1';env['OMP_NUM_THREADS']='1'
    env['PYTHONPATH']=str(root/'src')
    result=subprocess.run([sys.executable,'-c',code,str(omega),str(length),str(eta)],
        cwd=root,env=env,text=True,capture_output=True,check=True)
    z=json.loads(result.stdout.strip().splitlines()[-1])
    return complex(z['real'],z['imag'])


def pair(z):return {'real':float(np.real(z)),'imag':float(np.imag(z))}


def create(include_production=False):
    refs={};seq={}
    for w in FREQUENCIES:
        kr=keldysh_response(KSpec(omega=w,nodes_per_interval=640))
        refs[f'{w:.2f}']={'sigma_device':pair(kr['sigma_device']),
            'mixed_dB_sigma_device':pair(kr['mixed_dB_sigma_device'])}
        seq[f'{w:.2f}']=[]
        for L,eta in ((250,.02),(350,.02)):
            val=lehmann_response(Spec(omega=w,length=L,eta=eta))['sigma_device']
            seq[f'{w:.2f}'].append({'length':L,'eta':eta,'sigma_device':pair(val),
                'absolute_error':float(abs(val-kr['sigma_device']))})
    return {'version':'v70','scientific_status':'VERIFIED finite-L Lehmann/Ward identities; CONDITIONAL joint L/eta sequence; continuum discrepancy remains finite and original Xiang source remains PARTIAL',
        'model':'same device polarization-current observable; explicit finite tight-binding leads; independent full-Hamiltonian Lehmann response',
        'frequencies':list(FREQUENCIES),'keldysh_reference':refs,'quick_sequence':seq,
        'limit_path':{'medium':{'length':500,'eta':.01},'production':{'length':1000,'eta':.005},
            'warning':'L and eta must be varied jointly; intermediate errors need not decrease monotonically at every L due discrete finite-size resonances.',
            'eta_role':'eta is a spectral resolution for discrete transitions, not a phenomenological reservoir relaxation rate.'}}


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--output',default='results/v70_ordered_limits_quick.json')
    p.add_argument('--production-frequency',type=float,choices=(.20,.29,.40),
        help='run one L=1000, eta=.005 production point in an isolated invocation')
    a=p.parse_args();path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
    if a.production_frequency is None:
        data=create(False)
    else:
        w=a.production_frequency
        kr=keldysh_response(KSpec(omega=w,nodes_per_interval=640))
        value=production_value(w,length=1000,eta=.005)
        data={'version':'v70','frequency':w,'length':1000,'eta':.005,
            'sigma_device':pair(value),'keldysh_reference':pair(kr['sigma_device']),
            'absolute_error':float(abs(value-kr['sigma_device'])),
            'relative_error':float(abs(value-kr['sigma_device'])/abs(kr['sigma_device'])),
            'scientific_status':'finite-L/eta point; residual is nonzero and not claimed converged'}
    path.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print('v70 ordered-limit evidence written:',path)
