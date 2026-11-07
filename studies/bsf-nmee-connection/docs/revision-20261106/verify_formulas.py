"""Check the written Bloch blocks and interband formula against the existing model."""
from pathlib import Path
import sys,json
import numpy as np
repo=Path('/Users/jiaqi/Desktop/worksofJiaqi/quantum-geometry-reproductions')
sys.path[:0]=[str(repo/'studies/bsf-nmee-connection'),str(repo/'papers/jia-2026-nonlin-magnetoelectric-edelstein')]
from bsf_nmee.bilayer_audit import AuditBilayer,AuditSpec
from bsf_nmee.bilayer_protocol import ab_ba_shift
from nmee_lab.clean_core.models import I2,SX,SY,SZ
from nmee_lab.clean_core.core import chi_k
k=np.array([[.13,.29],[8.1,-3.2],[-5.8,12.2]])
a=.142;t0=.85;so=.010;rr=.020;m=.045;tp=.010;bias=.012
obj=AuditBilayer(AuditSpec(shift=ab_ba_shift(.25),staggered_eV=m,tunneling_eV=tp))
delta=a*np.array([[0,1],[np.sqrt(3)/2,-.5],[-np.sqrt(3)/2,-.5]])
bvec=np.array([delta[1]-delta[2],delta[2]-delta[0],delta[0]-delta[1]])
def mono(rashba):
 q=np.zeros((len(k),2,2),complex)
 for vec in delta:q-=np.exp(1j*k@vec)[:,None,None]*(t0*I2+1j*rashba*(SX*vec[1]-SY*vec[0])/a)
 dz=-2*so*np.sin(k@bvec.T).sum(1)
 h=np.zeros((len(k),4,4),complex)
 h[:,:2,:2]=m*I2+dz[:,None,None]*SZ;h[:,2:,2:]=-m*I2-dz[:,None,None]*SZ
 h[:,:2,2:]=q;h[:,2:,:2]=q.conj().swapaxes(1,2)
 return h
H=np.zeros((len(k),8,8),complex);H[:,:4,:4]=mono(rr)+bias/2*np.eye(4);H[:,4:,4:]=mono(-rr)-bias/2*np.eye(4)
H[:,:4,4:]=obj.tunneling(k);H[:,4:,:4]=obj.tunneling(k).conj().swapaxes(1,2)
ref=obj.matrices(k)+bias/2*np.kron(np.diag([1,-1]),np.eye(4));err=float(np.max(abs(H-ref)));assert err<1e-12
Hfield=H-.001*obj.spins['total_z'];e,u=np.linalg.eigh(Hfield);v=obj.matrices(k,axis=0);obs=np.broadcast_to(obj.spins['bottom_x'],H.shape)
occ=u[:,:,:4];emp=u[:,:,4:];w=(occ.conj().transpose(0,2,1)@obs@emp)*(emp.conj().transpose(0,2,1)@v@occ).transpose(0,2,1)
written=-2*(w/(e[:,4:][:,None,:]-e[:,:4,None])**2).imag.sum((1,2))
reference=chi_k(e,u,v,obs,np.tile([1,1,1,1,0,0,0,0],(len(k),1)))
err2=float(np.max(abs(written-reference)));assert err2<1e-12
L=np.kron(np.diag([1,-1]),np.eye(4));neg=H.copy();neg[:,:4,4:]*=-1;neg[:,4:,:4]*=-1
err3=float(np.max(abs(neg-L@H@L)));assert err3<1e-12
result={'Bloch_blocks_max_eV':err,'interband_formula_max':err2,'t_sign_unitary_max_eV':err3,'scope':'Algebra at three off-symmetry k points; no new integrated response or DFT.'}
Path(__file__).with_name('formula_checks.json').write_text(json.dumps(result,indent=2));print(result)
