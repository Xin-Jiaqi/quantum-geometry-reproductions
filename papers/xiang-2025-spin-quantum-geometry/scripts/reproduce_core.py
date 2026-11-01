"""Four-band full spin-current response, independent of reduced Eq13."""
import sys,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from methods.response import SpinModel,integrate_response

def run():
    rng=np.random.default_rng(73)
    def matrix():
        a=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4));return .08*(a+a.conj().T)
    A,B,C=matrix(),matrix(),matrix()
    h0=np.diag([-1.3,-.3,.5,1.4])+matrix()
    D=matrix()  # Independent sin(kx) term breaks the shifted inversion of the first three harmonics.
    sz=np.kron(np.eye(2),np.diag([1.,-1.]))
    model=SpinModel(lambda x,y:h0+np.cos(x)*A+np.sin(y)*B+np.sin(x+y)*C+np.sin(x)*D,
                    lambda x,y,a:np.zeros((4,4)),lambda x,y,c:sz,2*np.pi*np.eye(2))
    settings=dict(energy1=.08+.12j,energy2=-.08+.12j,axis1=0,axis2=1,
                  mu=.1,thermal_energy=.25,current_axis=0,component='z',step=.003)
    rows=[]
    for n,h in ((25,.003),(41,.003),(41,.006)):
        q=-np.pi+2*np.pi*(np.arange(n)+.31)/n
        points=np.stack(np.meshgrid(q,q,indexing='ij'),-1).reshape(-1,2)
        r=integrate_response(model,points,**dict(settings,step=h))
        rows.append(dict(n=n,step=h,response=r['response']))
    mesh_error=abs(rows[1]['response']-rows[0]['response'])
    step_error=abs(rows[1]['response']-rows[2]['response'])
    if abs(rows[1]['response'])<1e-6:raise AssertionError('Need a nonzero transport benchmark')
    if max(mesh_error,step_error)>2e-5:
        raise AssertionError(f'Convergence target not reached: mesh={mesh_error}, step={step_error}, rows={rows}')
    return dict(protocol='full_quadratic_spin_current_symmetric',model='synthetic four-band physical spin',
                settings=settings,rows=rows,mesh_change=mesh_error,step_change=step_error,
                units='e=hbar=1; spin sigma/2; d2k/(2pi)^2',
                scope='Not reduced Eq13, not isolated tau0, not a DFT material')

if __name__=='__main__':print(json.dumps(run(),indent=2,default=lambda z:dict(real=z.real,imag=z.imag)))
