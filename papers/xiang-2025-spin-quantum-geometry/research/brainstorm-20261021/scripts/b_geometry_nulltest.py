"""Minimal Berry-curvature test of pure layer rephasing versus modulated hopping.
Uses occupied two-band P=(1-n.sigma)/2 and F=i Tr P[d_k P,d_delta P].
No claim that these toy Hamiltonians represent a particular real material.
"""
import numpy as np
nk=3000
k=np.linspace(-np.pi,np.pi,nk,endpoint=False)

def curvature(d,dk,dd):
    return np.einsum('ij,ij->i',d,np.cross(dk,dd))/(2*np.linalg.norm(d,axis=1)**3)

# 1) A parameter that merely rephases layer orbitals: H(delta)=exp(-i delta tau_z/2)H(0)exp(i delta tau_z/2).
dx=0.7+0.2*np.cos(k); dy=0.25*np.sin(k); dz=0.5+0.3*np.cos(k)
d=np.column_stack((dx,dy,dz))
dk=np.column_stack((-0.2*np.sin(k),0.25*np.cos(k),-0.3*np.sin(k)))
ddel=np.column_stack((-dy,dx,np.zeros_like(k)))
F=curvature(d,dk,ddel)
nhat_z=dz/np.linalg.norm(d,axis=1)
F_id=-np.gradient(nhat_z,2*np.pi/nk,edge_order=2)/2
print('REPHASING gauge-only k-integral F:',F.mean())
print('REPHASING max local |F|:',np.max(np.abs(F)))
print('REPHASING max analytic derivative error (bulk):',np.max(np.abs(F[2:-2]-F_id[2:-2])))

# 2) A structural parameter changes actual hopping magnitude and onsite imbalance.
# Rice--Mele-like toy: d=(t1(delta)+t2 cos k, t2 sin k, m sin delta).
# δ is a generic cyclic structural coordinate, not claimed a realistic sliding path.
t0=1.10; amp=0.50; t2=0.90; mass=0.35
for delta in [0.2,0.7,1.3,2.0,3.5,5.4]:
    t1=t0+amp*np.cos(delta)
    dr=np.column_stack((t1+t2*np.cos(k),t2*np.sin(k),mass*np.sin(delta)*np.ones_like(k)))
    drk=np.column_stack((-t2*np.sin(k),t2*np.cos(k),np.zeros_like(k)))
    drd=np.column_stack((-amp*np.sin(delta)*np.ones_like(k),np.zeros_like(k),mass*np.cos(delta)*np.ones_like(k)))
    fm=curvature(dr,drk,drd)
    print('REAL MOD delta %.2f  k-integral F %.9f  gap %.5f'%(delta,fm.mean(),2*np.min(np.linalg.norm(dr,axis=1))))
# Integrate mixed curvature over torus (k,delta) to find Chern number.
Ndelta=801
vals=[]
for delta in np.linspace(-np.pi,np.pi,Ndelta,endpoint=False):
    t1=t0+amp*np.cos(delta)
    dr=np.column_stack((t1+t2*np.cos(k),t2*np.sin(k),mass*np.sin(delta)*np.ones_like(k)))
    drk=np.column_stack((-t2*np.sin(k),t2*np.cos(k),np.zeros_like(k)))
    drd=np.column_stack((-amp*np.sin(delta)*np.ones_like(k),np.zeros_like(k),mass*np.cos(delta)*np.ones_like(k)))
    vals.append(curvature(dr,drk,drd).mean())
print('REAL MOD k-delta Chern:',2*np.pi*np.mean(vals))
