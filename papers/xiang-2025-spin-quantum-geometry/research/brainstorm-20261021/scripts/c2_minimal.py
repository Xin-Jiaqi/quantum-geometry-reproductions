import numpy as np
from itertools import product
I=np.eye(2, dtype=complex);sx=np.array([[0,1],[1,0]],complex);sy=np.array([[0,-1j],[1j,0]],complex);sz=np.diag([1,-1]).astype(complex)
def X(a,b):return np.kron(a,b)
A=X(sx,sx);PT=X(sx,1j*sy);C2=X(I,1j*sz);S=X(sz,sx)/2
# H = t tau_x + u cos kx tau_y + lam tau_y sigma_x + d sin kx tau_z + M tau_z sigma_z + beta sin ky tau_z sigma_x
p=dict(t=1.0,u=.21,lam=.33,d=.49,M=1.3,beta=.24)

def ham(x,y,q=1.0):
 t,u,lam,d,M,beta=[p[z] for z in ('t','u','lam','d','M','beta')]
 hp=q*(t*X(sx,I)+u*np.cos(x)*X(sy,I)+lam*X(sy,sx))
 hlocal=d*np.sin(x)*X(sz,I)+M*X(sz,sz)+beta*np.sin(y)*X(sz,sx)
 return hp+hlocal,hp

def vel(x,y,q=1.0):
 return -q*p['u']*np.sin(x)*X(sy,I)+p['d']*np.cos(x)*X(sz,I)

def k_resp(N,w,q=1.,fd=False):
 ps=np.linspace(-np.pi,np.pi,N,endpoint=False);acc=np.zeros(3,complex);gap=100
 for kx,ky in product(ps,ps):
  h,hp=ham(kx,ky,q);en,U=np.linalg.eigh(h);gap=min(gap,en[2]-en[1]);
  O=[1j*(hp@S-S@hp), 1j*((h-hp)@S-S@(h-hp)), S]
  if fd:
   eps=1.e-5
   def P(xx):
    _,v=np.linalg.eigh(ham(xx,ky,q)[0]);return v[:,:2]@v[:,:2].conj().T
   dP=(P(kx+eps)-P(kx-eps))/(2*eps)
   cov=1j*(U.conj().T@dP@U)
   R=np.zeros((4,4),complex)
   for i in range(4):
    for j in range(4):
     if (i<2)!=(j<2):R[i,j]=cov[i,j]/(w-en[i]+en[j])
   for a,observable in enumerate(O):acc[a]+=np.trace((U.conj().T@observable@U)@R)/N**2
  else:
   v=U.conj().T@vel(kx,ky,q)@U
   for j,observable in enumerate(O):
    oo=U.conj().T@observable@U
    for a in range(2):
     for b in range(2,4):
      D=en[b]-en[a];Q=oo[a,b]*(-1j*v[b,a]/D)
      acc[j]+=(Q/(w-D)-Q.conj()/(w+D))/N**2
 return gap,acc
x,y=.234,.742;h,hp=ham(x,y)
print('A symmetry',np.max(np.abs(A@h.conj()@A.conj().T-ham(-x,-y)[0])))
print('PT violation',np.max(np.abs(PT@h.conj()@PT.conj().T-h)))
print('C2 violation',np.max(np.abs(C2@h@C2.conj().T-ham(-x,-y)[0])))
for w in (0,.1,.2,.3):
 gap,r=k_resp(39,w)
 print('w',w,'gap',gap,'Tint/Tloc/S',r,'Ward',r[0]+r[1], 'ratio_to_-iws', (r[0]+r[1])/(-1j*w*r[2]) if w else '-')
print('fd check',k_resp(25,.2),k_resp(25,.2,fd=True))
for q in (1,.5,.2,.0):
 gap,r=k_resp(29,.2,q)
 print('hopping scale',q,'gap',gap,'Tint',r[0],'S',r[2])

def curvature_mixed(N=21, hstep=1.e-4,theta_step=1.e-4):
 kpts=np.linspace(-np.pi,np.pi,N,endpoint=False)
 result=0.
 # finite difference P(k,theta), perturbing physical interlayer hoppings only
 ee,UU=np.linalg.eigh(S)
 def spinrot(z):return UU@np.diag(np.exp(-1j*z*ee))@UU.conj().T
 def proj(kx,ky,theta=0):
  hh,hp=ham(kx,ky)
  rot=spinrot(theta)
  newH=hh-hp+rot@hp@rot.conj().T
  _,U=np.linalg.eigh(newH)
  return U[:,:2]@U[:,:2].conj().T
 for x,y in product(kpts,kpts):
  p0=proj(x,y)
  dk=(proj(x+hstep,y)-proj(x-hstep,y))/(2*hstep)
  dt=(proj(x,y,theta_step)-proj(x,y,-theta_step))/(2*theta_step)
  result+= (1j*np.trace(p0@(dk@dt-dt@dk))).real/N**2
 return result
if __name__=='__main__':
 print('mixed curvature integral',curvature_mixed(21), 'static Kubo',k_resp(21,0)[1][0])