"""v53: authentic orbital SM Eq(49)->(52) metric and three-band, Eq(56) PHE.

Primary source the user's 16-page Supplemental Material, pp.11-12.
Eq49/F has explicit factor 2 in band-velocity*metric term, while printed
Eq52 has factor 1. Test with a true periodic 3D two-level Hamiltonian
whose local interband metric is 1/4 d_i(nhat).d_j(nhat).
Compute the moment-curvature term in a genuine periodic 3-band lattice
and separate m=l from m!=l. Eq56's first/second lines are checked
independently by periodic integration by parts in a two-level Bloch model.
No physical transport/bath completion or author manuscript correction implied.
"""
from __future__ import annotations
import numpy as np
SM_PDF_SHA256='9b0b8603d13962e100cc059bf4d96b91aa37062360d13decd1cb2eda495332c4'
MAIN_PDF_SHA256='b87e0970f25dd45c5048f8cbbb8691b53ea215040935e9b21a4cde80c69b2730'

def _check_N(N):
    if not isinstance(N,int) or N<12: raise ValueError('periodic mesh N>=12')

def metric_eq49_to_eq52(N=32,mass=2.4,mu=2.55,T=.34,tilt_scale=1.):
  _check_N(N)
  k=-np.pi+(np.arange(N)+.39)*2*np.pi/N
  x,y,z=np.meshgrid(k,k,k,indexing='ij')
  d0=tilt_scale*(.3*np.sin(x)+.27*np.sin(y)+.2*np.sin(z)+.11*np.cos(y))
  d0gr=tilt_scale*np.array([.3*np.cos(x),.27*np.cos(y)-.11*np.sin(y),.2*np.cos(z)])
  dx=.55+np.sin(x)+.3*np.cos(y)+.18*np.sin(z)
  dy=.23+np.sin(y)+.31*np.sin(z)+.12*np.cos(x)
  dz=mass+.42*np.cos(x)+.23*np.cos(y)+.37*np.cos(z)+.25*np.sin(x)
  d=np.array([dx,dy,dz]);D=np.linalg.norm(d,axis=0);u=d/D
  dg=np.array([
     [np.cos(x),-.12*np.sin(x),.25*np.cos(x)-.42*np.sin(x)],
     [-.3*np.sin(y),np.cos(y),-.23*np.sin(y)],
     [.18*np.cos(z),.31*np.cos(z),-.37*np.sin(z)]
  ]) # (deriv i, orbital component)
  du=np.zeros_like(dg)
  for i in range(3):
    alpha=np.sum(u*dg[i],axis=0)
    du[i]=(dg[i]-u*alpha)/D
  g=np.einsum('icxyz,jcxyz->ijxyz',du,du)*.25
  h=2*np.pi/N
  def der(a,axis): return (np.roll(a,-1,axis=axis)-np.roll(a,1,axis=axis))/(2*h)
  # eps_yzx=+; eps_yxz=-, eps_xyz=+; eps_xzy=-
  grad_g_xz=np.array([der(g[0,2],i) for i in range(3)])
  grad_g_yz=np.array([der(g[1,2],i) for i in range(3)])
  grad_g_zz=np.array([der(g[2,2],i) for i in range(3)])
  accum=np.zeros(4)
  for sign in (-1,1):
    E=d0+sign*D;gap=2*sign*D
    v0=d0gr
    dD=np.einsum('cxyz,icxyz->ixyz',u,dg)
    vn=d0gr+sign*dD
    vsum=2*d0gr
    f=1/(1+np.exp(np.clip((E-mu)/T,-100,100)))
    fp=-f*(1-f)/T
    A_zx=2*vsum[2]*g[0,2]/gap + grad_g_xz[2]
    A_xz=2*vsum[0]*g[2,2]/gap + grad_g_zz[0]
    A_yz=2*vsum[1]*g[2,2]/gap + grad_g_zz[1]
    A_zy=2*vsum[2]*g[1,2]/gap + grad_g_yz[2]
    F_yz=-A_zx+A_xz
    F_xz=-A_yz+A_zy
    eq48=(der(F_yz,0)-der(F_xz,1))*f
    vterm=(vn[0]*(vsum[2]*g[0,2]-vsum[0]*g[2,2])/gap-vn[1]*(vsum[1]*g[2,2]-vsum[2]*g[1,2])/gap)*fp
    hess=(der(grad_g_xz[2],0)-der(grad_g_zz[0],0)-der(grad_g_zz[1],1)+der(grad_g_yz[2],1))*f
    eq52=vterm-hess
    eq52_correct=2*vterm-hess
    accum+=np.array([eq48.mean(),eq52.mean(),eq52_correct.mean(),vterm.mean()])
  return {'N':N,'source48_F':accum[0],'printed52_metric':accum[1],'factor2_corrected':accum[2],'missing_velocity':accum[3], 'residual_original':accum[0]-accum[1], 'residual_correct':accum[0]-accum[2]}

def orbital_moment_three_band(N=56,mu=.5,T=.5,band_phases=(0.,0.,0.)):
  _check_N(N)
  if len(band_phases)!=3 or not np.isfinite(band_phases).all(): raise ValueError("three finite U1 band phases")
  k=-np.pi+(np.arange(N)+.33)*2*np.pi/N
  x,y=np.meshgrid(k,k,indexing='ij')
  H=np.zeros((N,N,3,3),complex)
  H[...,0,0]=-1.3+.23*np.cos(x)
  H[...,1,1]=.3+.3*np.cos(y)
  H[...,2,2]=2.1+.2*np.cos(x+y)
  h01=.43*(np.sin(x)+1j*np.sin(y))+.13j*np.cos(x)
  h02=.38*(np.sin(y)+1j*np.sin(x))
  h12=.49*(np.sin(x+y)+1j*np.cos(x))
  for (i,j),h in [((0,1),h01),((0,2),h02),((1,2),h12)]:
     H[...,i,j]=h;H[...,j,i]=h.conj()
  vx=np.zeros_like(H);vy=np.zeros_like(H)
  vx[...,0,0]=-.23*np.sin(x);vy[...,1,1]=-.3*np.sin(y)
  vx[...,2,2]=-.2*np.sin(x+y);vy[...,2,2]=-.2*np.sin(x+y)
  pairs=[(0,1,.43*np.cos(x)-.13j*np.sin(x),1j*.43*np.cos(y)),(0,2,1j*.38*np.cos(x),.38*np.cos(y)),(1,2,.49*(np.cos(x+y)-1j*np.sin(x)),.49*np.cos(x+y))]
  for i,j,a,b in pairs:
   vx[...,i,j]=a;vx[...,j,i]=a.conj()
   vy[...,i,j]=b;vy[...,j,i]=b.conj()
  e,u=np.linalg.eigh(H)
  u=u*np.exp(1j*np.asarray(band_phases))[None,None,None,:]
  ud=u.conj().swapaxes(-1,-2)
  vxb=ud@vx@u;vyb=ud@vy@u
  om=np.zeros((N,N,3,3),float)
  for n in range(3):
   for m in range(3):
    if n==m:continue
    denom=(e[...,m]-e[...,n])**2
    # r_nm = i v_nm/(Em-En), Ωnm = i(rnx rmy-rny rmx)= -2 Im(rnx*rmy)
    om[...,n,m]=np.real(1j*(vxb[...,n,m]*vyb[...,m,n]-vyb[...,n,m]*vxb[...,m,n])/denom)
  fp=-1/(4*T*np.cosh(np.clip((e-mu)/(2*T),-50,50))**2)
  full=0;diag=0
  per=[]
  for n in range(3):
   eps=e[...,n,None]-e
   allomega=om[...,n,:].sum(axis=-1)
   mw=(eps*om[...,n,:]).sum(axis=-1)
   fullterm=-.5*np.mean(fp[...,n]*allomega*mw)
   diagterm=-.5*np.mean(fp[...,n]*(eps*om[...,n,:]**2).sum(axis=-1))
   per.append((fullterm,diagterm,fullterm-diagterm))
   full+=fullterm;diag+=diagterm
  return {'N':N,'mu':mu,'T':T,'complete_moment_Berry':full,'printed_eq52_m_eq_l':diag,'omitted_threeband':full-diag,'ratio_abs_to_complete':abs((full-diag)/full) if full else None,'per_band':per,'max_om':np.max(abs(om))}

def orbital_phe_eq56(N=32,mu=2.55,T=.36):
    _check_N(N)
    k=-np.pi+(np.arange(N)+.39)*2*np.pi/N
    x,y=np.meshgrid(k,k,indexing='ij')
    d0=.36*np.sin(x)+.21*np.sin(y)+.12*np.cos(x+y)
    dx=.4+np.sin(x)+.23*np.cos(y)
    dy=.14+np.sin(y)+.37*np.cos(x)
    dz=2.5+.34*np.cos(x)+.25*np.cos(y)+.17*np.sin(x+y)
    d=np.array([dx,dy,dz]);R=np.linalg.norm(d,axis=0);u=d/R
    dgx=np.array([np.cos(x),-.37*np.sin(x),-.34*np.sin(x)+.17*np.cos(x+y)])
    dgy=np.array([-.23*np.sin(y),np.cos(y),-.25*np.sin(y)+.17*np.cos(x+y)])
    rx=np.sum(u*dgx,axis=0);ry=np.sum(u*dgy,axis=0)
    pref=np.sum(u*np.cross(dgx,dgy,axis=0),axis=0)/(2*R*R)
    h=2*np.pi/N
    deriv=lambda a,i:(np.roll(a,-1,axis=i)-np.roll(a,1,axis=i))/(2*h)
    total1=total2=orb=0
    for s in (-1,1):
        E=d0+s*R
        vx=.36*np.cos(x)-.12*np.sin(x+y)+s*rx
        vy=.21*np.cos(y)-.12*np.sin(x+y)+s*ry
        Om=s*pref
        m=R*pref
        f=1/(1+np.exp(np.clip((E-mu)/T,-100,100)))
        fp=-f*(1-f)/T
        fpp=fp*(-(1-2*f)/T)
        hess_f=fpp*vx*vy+fp*deriv(vx,1)
        termA=(-deriv(m,0)*vy*fp-vx*deriv(m*fp,1)).mean()
        termB=(2*m*hess_f-vx*vy*m*fpp).mean()
        curvature=(-vx*vy*Om*fp).mean()
        total1+=termA+curvature
        total2+=termB+curvature
        orb+=curvature
    return {'N':N,'first_eq56':total1,'second_eq56':total2,'absolute_delta':total1-total2,'curvature_term':orb}

def full_audit():
    return {
        'scientific_status':'PARTIAL',
        'primary_source':'original SM pages 11–12, Eqs. (48),(49),(52),(56)',
        'original_source_SHA256':{'sm':SM_PDF_SHA256,'main':MAIN_PDF_SHA256},
        'metric_eq49_52_mesh':[metric_eq49_to_eq52(N) for N in (20,28,36,48)],
        'three_band_curvature_moment':[orbital_moment_three_band(N,mu=-1.1,T=.45) for N in (32,56)],
        'orbital_phe_eq56_mesh':[orbital_phe_eq56(N) for N in (24,32,48,64)],
        'scope':{
          'eq49_52':'The metric-velocity term printed with coefficient 1 in Eq52 does not equal Eq49 after integration by parts; coefficient 2 is supported in a two-band physical Bloch metric.',
          'threeband':'Eq52 explicitly drops m!=l of the final Berry-curvature–orbital-moment sector; this term can be order-one in an independent three-band model. Not a total-MNHE percentage.',
          'eq56':'Two printed representations are consistent after BZ integration; residual decreases as finite difference mesh is refined.',
          'limits':['Original Eq49 already discards a separate class of 3-band terms',
             'No claim of complete experimental or conserving orbital response',
             'SM Eq20/31 sign and Fig2 SI e2/h vs e2/hbar remain open']}}
