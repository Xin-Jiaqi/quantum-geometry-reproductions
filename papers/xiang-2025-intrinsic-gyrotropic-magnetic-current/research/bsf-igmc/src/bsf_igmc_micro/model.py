"""A distinct eight-band bilayer real-space hopping study, not Xiang's model.
Basis: (top,bottom) x (orbital A,B) x (spin up,down).
Onsite electrostatic layer-odd potentials are EXCLUDED by construction.
The two orbitals are at different in-plane y coordinates; all interlayer
hoppings are Gaussian overlaps between explicitly displaced Wannier centers.
"""
from dataclasses import dataclass
import numpy as np
I=np.eye(2,dtype=complex)
X=np.array([[0,1],[1,0]],complex)
Y=np.array([[0,-1j],[1j,0]],complex)
Z=np.diag([1.,-1.]).astype(complex)
kron=lambda *a:__import__('functools').reduce(np.kron,a)
S={a:kron(I,I,b) for a,b in [('x',X),('y',Y),('z',Z)]}
LZ=kron(Z,I,I)
MZ=kron(X,I,Z)
MZT=MZ@kron(I,I,1j*Y)
MXT=kron(I,I,Z)

@dataclass(frozen=True)
class Parameters:
  yB:float=.37
  xB:float=0.
  t0:float=.45
  kappa:float=2.8
  tAB:float=.22
  mass:float=.42
  lambda_so:float=.38
  lambda_z:float=.15
  exchange:float=.23
  cant:float=.14
  inter_soc:float=.16
  inter_spin:float=.18
  shell:int=2
  inter_orbital:float=1.

def overlaps(kx,ky,u,p):
  """t_{ab}(R,u) ~ exp[-κ|R+r_b-r_a+(0,u)|²].
  Return T, exact ∂kx T and ∂ky T; no fitted odd-u onsite terms.
  spin-independent, SOC, and spin-z hopping use same geometric overlap.
  """
  yy=np.array([0.,p.yB]);xx=np.array([0.,p.xB])
  t=np.zeros((2,2),complex);dx=t.copy();dy=t.copy()
  for rx in range(-p.shell,p.shell+1):
    for ry in range(-p.shell,p.shell+1):
      phase=np.exp(1j*(kx*rx+ky*ry))
      distx=rx+xx[None,:]-xx[:,None]
      disty=ry+yy[None,:]-yy[:,None]+u
      amp=p.t0*np.exp(-p.kappa*(distx**2+disty**2))*phase
      weight=np.array([[1.,p.inter_orbital],[p.inter_orbital,1.]])
      amp*=weight
      t+=amp;dx+=1j*rx*amp;dy+=1j*ry*amp
  return t,dx,dy

def hamiltonian(kx,ky,u,branch='AF',p=Parameters(),onsite_layer=0.):
  if branch not in ('AF','FM'):raise ValueError(branch)
  if onsite_layer!=0:raise ValueError('Independent layer-odd onsite forbidden: use geometric hopping only')
  # Orbital onsite/hybridization is identical on both layers; fixed material chemistry.
  orb=p.mass*Z+2*p.tAB*np.cos(ky)*X
  d_orb_y=-2*p.tAB*np.sin(ky)*X
  # Spin-dependent real-space hopping is a two-center SOC ansatz; all its
  # (k,u) dependence derives from the same geometrical overlap function.
  t,tx,ty=overlaps(kx,ky,u,p)
  def block(v):
    return kron(v,I)-1j*p.inter_soc*kron(v,X)+(p.inter_spin if branch=='AF' else 0)*kron(v,Z)
  upper=block(t)
  H=np.zeros((8,8),complex)
  H[:4,4:]=upper;H[4:,:4]=upper.conj().T
  deriv=[]
  for v in (tx,ty):
    d=np.zeros((8,8),complex);d[:4,4:]=block(v);d[4:,:4]=block(v).conj().T;deriv.append(d)
  base=kron(I,orb,I)
  rashba=p.lambda_so*kron(Z,I,np.sin(ky)*X-np.sin(kx)*Y)
  magnetic=p.exchange*kron(Z if branch=='AF' else I,I,Y)
  cant=p.cant*kron(I if branch=='AF' else Z,I,Z)
  intra=p.lambda_z*np.sin(kx)*kron(I,I,Z)
  H+=base+rashba+magnetic+cant+intra
  vx=deriv[0]-p.lambda_so*np.cos(kx)*kron(Z,I,Y)+p.lambda_z*np.cos(kx)*kron(I,I,Z)
  vy=deriv[1]+kron(I,d_orb_y,I)+p.lambda_so*np.cos(ky)*kron(Z,I,X)
  return H,{'x':vx,'y':vy}

def connection_error(u,branch,p=Parameters(),N=5):
  result={'partner':0.,'internal_mxT':0.,'hermiticity':0.,'derivative_x':0.,'derivative_y':0.}
  for kx in np.linspace(-2,2,N):
    for ky in np.linspace(-2.1,2.2,N):
      H,V=hamiltonian(kx,ky,u,branch,p)
      if branch=='AF':
        Q,_=hamiltonian(kx,ky,-u,branch,p)
        expected=MZ@H@MZ.conj().T
      else:
        Q,_=hamiltonian(-kx,-ky,-u,branch,p)
        expected=MZT@H.conj()@MZT.conj().T
      I_H,_=hamiltonian(kx,-ky,u,branch,p)
      result['partner']=max(result['partner'],float(np.max(np.abs(Q-expected))))
      result['internal_mxT']=max(result['internal_mxT'],float(np.max(np.abs(I_H-MXT@H.conj()@MXT.conj().T))))
      result['hermiticity']=max(result['hermiticity'],float(np.max(np.abs(H-H.conj().T))))
      eps=2e-5
      for a,q in [('x',0),('y',1)]:
        kp=[kx,ky];km=kp.copy();kp[q]+=eps;km[q]-=eps
        Hp,_=hamiltonian(*kp,u,branch,p);Hm,_=hamiltonian(*km,u,branch,p)
        err=float(np.max(np.abs((Hp-Hm)/(2*eps)-V[a])))
        result['derivative_'+a]=max(result['derivative_'+a],err)
  return result

def geometric_response(u,branch='AF',p=Parameters(),N=31,alpha='x',beta='x',mu=0.,strict_gap=True):
  """Bloch eigenbasis D_xx/ω = Σ_{occ,empty} -2 Im(v_nm s_mn)/Δ².
  All bands and identical same-gauge spin operators included.
  The absolute normalization is the explicit model-only BZ average.
  """
  ks=-np.pi+(np.arange(N)+.5)*2*np.pi/N
  accum=0.;ecmax=-1e10;evmin=1e10;below=-1e10;above=1e10
  occupancy_min=8;occupancy_max=0;fs_c=0.
  for kx in ks:
    for ky in ks:
      H,v=hamiltonian(kx,ky,u,branch,p)
      E,U=np.linalg.eigh(H)
      occ=E<mu;cnt=int(np.sum(occ));occupancy_min=min(cnt,occupancy_min);occupancy_max=max(cnt,occupancy_max)
      # At half filling require exactly 4 occupied bands and global indirect gap.
      below=max(below,float(E[3])); above=min(above,float(E[4]))
      A=U.conj().T@v[alpha]@U;B=U.conj().T@S[beta]@U
      for n in range(8):
        for m in range(n+1,8):
          f=float(occ[n])-float(occ[m])
          if f==0:continue
          gap=E[m]-E[n];prod=A[n,m]*B[m,n]
          accum+=f*(-2*np.imag(prod)/gap**2)
          fs_c+=f*(-2*np.real(prod)/gap)
  if strict_gap and (occupancy_min!=4 or occupancy_max!=4 or above-below<1e-5):
    raise ValueError(f'not an indirect half-filled insulator at μ={mu}: occupancy={occupancy_min}-{occupancy_max}, gap={above-below}')
  return {'Dxx_over_omega':float(accum/N**2),'Cxx':float(fs_c/N**2),
          'indirect_gap':float(above-below),'E4_max':below,'E5_min':above,'occ_range':[occupancy_min,occupancy_max],'mesh':N}

def layer_polarization(u,branch='AF',p=Parameters(),N=31,mu=0.):
  ks=-np.pi+(np.arange(N)+.5)*2*np.pi/N
  accum=0.
  for kx in ks:
    for ky in ks:
      H,_=hamiltonian(kx,ky,u,branch,p)
      E,U=np.linalg.eigh(H)
      occ=U[:,E<mu]
      accum+=float(np.trace(occ.conj().T@LZ@occ).real)
  return accum/N**2

def wilson_zeeman_slope(u,branch='AF',p=Parameters(),N=31,field=5e-4):
  """Independent Wilson determinant dP_x/dB_x from four occupied bands.
  Periodic-gauge Berry phase in x; orbital x positions are identical (=0).
  Magnetic perturbation H(B)=H-B*s_x; compares to Dxx/omega kernel.
  """
  ks=-np.pi+(np.arange(N)+.5)*2*np.pi/N
  phi=[]
  for bz in (-field,field):
    eigs=np.empty((N,N,8,4),complex)
    for i,kx in enumerate(ks):
      for j,ky in enumerate(ks):
        H,_=hamiltonian(kx,ky,u,branch,p)
        e,U=np.linalg.eigh(H-bz*S['x'])
        if e[3]>=0 or e[4]<=0:raise ValueError('Zeeman perturbation gap closed')
        eigs[i,j]=U[:,:4]
    overlaps=np.einsum('...ia,...ib->...ab',eigs.conj(),np.roll(eigs,-1,axis=0))
    determinants=np.linalg.det(overlaps)
    if np.min(np.abs(determinants))<1e-8:raise ValueError('singular Wilson link')
    phi.append(-np.angle(np.prod(determinants/np.abs(determinants),axis=0)))
  dphi=np.angle(np.exp(1j*(phi[1]-phi[0])))
  return float(np.sum(dphi)/(2*field*N*2*np.pi))

def band_energy(u,branch='AF',p=Parameters(),N=31):
  """Occupied band-energy SUM only. NOT total energy, no ionic repulsion!
  """
  ks=-np.pi+(np.arange(N)+.5)*2*np.pi/N
  accum=0.
  for kx in ks:
    for ky in ks:
      H,_=hamiltonian(kx,ky,u,branch,p)
      e=np.linalg.eigvalsh(H)
      accum+=float(e[:4].sum())
  return float(accum/N**2)
