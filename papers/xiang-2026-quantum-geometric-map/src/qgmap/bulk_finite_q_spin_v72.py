"""v72: periodic finite-q scalar-potential/spin-Zeeman mixed response.
The longitudinal current is exact [H(k+q)-H(k)]/q; transverse current a
DECLARED midpoint vertex. Spin source B*sigma_x, not orbital Peierls.
This does NOT reconstruct source selected 1/eta_e pole or original dc PHE.
"""
from dataclasses import dataclass,replace
import numpy as np

I=np.eye(2,dtype=complex)
SX=np.array([[0,1],[1,0]],complex)
SY=np.array([[0,-1j],[1j,0]],complex)
SZ=np.diag([1.,-1.]).astype(complex)

@dataclass(frozen=True)
class Spec:
    n:int=24
    q_steps:int=1
    omega:float=.27
    eta:float=.045
    field:float=.11
    mu:float=.32
    temperature:float=.12
    beta:float=.24
    alpha:float=.13
    mass:float=.32
    wilson:float=.35
    orbital_gauge:float=0.
    def __post_init__(self):
        vals=[self.omega,self.eta,self.field,self.mu,self.temperature,self.beta,
              self.alpha,self.mass,self.wilson,self.orbital_gauge]
        if not np.isfinite(vals).all() or self.eta<=0 or self.temperature<=0:
            raise ValueError('finite parameters and eta,T > 0 required')
        if not isinstance(self.n,int) or self.n<8 or not isinstance(self.q_steps,int) or not 1<=self.q_steps<self.n//2:
            raise ValueError('n>=8 and integer 1 <= q_steps < n//2')

def fermi(e,mu,T):
    x=np.clip((np.asarray(e)-mu)/T,-700,700)
    return 1/(1+np.exp(x))
def _dag(a):return np.swapaxes(a.conj(),-1,-2)
def _tr(a):return np.trace(a,axis1=-2,axis2=-1)

def _field_and_vertices(s):
    k=2*np.pi*np.arange(s.n)/s.n
    x,y=np.meshgrid(k,k,indexing='ij')
    w=2-np.cos(x)-np.cos(y)
    h=((s.beta*np.sin(x)+s.alpha*w)[...,None,None]*I
       -np.sin(y)[...,None,None]*SX+np.sin(x)[...,None,None]*SY
       +(s.mass+s.wilson*w)[...,None,None]*SZ+s.field*SX)
    vx=((s.beta*np.cos(x)+s.alpha*np.sin(x))[...,None,None]*I
        +np.cos(x)[...,None,None]*SY+(s.wilson*np.sin(x))[...,None,None]*SZ)
    z=np.broadcast_to(SX,h.shape).copy()
    if s.orbital_gauge:
        theta=s.orbital_gauge*(np.sin(x)+.3*np.cos(y))
        U=np.zeros_like(h)
        U[...,0,0]=np.exp(1j*theta)
        U[...,1,1]=np.exp(-1j*theta)
    else:U=np.broadcast_to(I,h.shape).copy()
    return h,vx,z,U

def _entries(s):
    h,vx,z,U=_field_and_vertices(s)
    shift=lambda x:np.roll(x,-s.q_steps,axis=1)
    hp,hm=shift(h),h
    vp,vm=shift(vx),vx
    zp,zm=shift(z),z
    up,um=shift(U),U
    hp_t=_dag(up)@hp@up
    hm_t=_dag(um)@hm@um
    zp_t=_dag(up)@zp@up
    zm_t=_dag(um)@zm@um
    W=_dag(up)@um
    Jx=_dag(up)@(vp+vm)@um*.5
    q=2*np.pi*s.q_steps/s.n
    Jy=_dag(up)@(hp-hm)@um/q
    return hp_t,hm_t,zp_t,zm_t,W,Jx,Jy,q

def _frechet_f(e,u,z,s):
    f=fermi(e,s.mu,s.temperature)
    d=e[..., :,None]-e[...,None,:]
    fd=f[..., :,None]-f[...,None,:]
    mid=fermi((e[..., :,None]+e[...,None,:])/2,s.mu,s.temperature)
    slope=-mid*(1-mid)/s.temperature
    mask=abs(d)>1e-10
    dd=np.divide(fd,d,out=np.zeros_like(fd),where=mask)
    dd=np.where(mask,dd,slope)
    return u@(dd*(_dag(u)@z@u))@_dag(u)

def _matrix_data(s):
    hp,hm,zp,zm,W,Jx,Jy,q=_entries(s)
    ep,up=np.linalg.eigh(hp);em,um=np.linalg.eigh(hm)
    fp=fermi(ep,s.mu,s.temperature);fm=fermi(em,s.mu,s.temperature)
    Fp=(up*fp[...,None,:])@_dag(up)
    Fm=(um*fm[...,None,:])@_dag(um)
    src=W@Fm-Fp@W
    sb=_dag(up)@src@um
    de=ep[..., :,None]-em[...,None,:]
    den=s.omega+1j*s.eta-de
    rb=sb/den
    R=up@rb@_dag(um)
    return locals()

def _response(R,W,Jx,Jy):
    v=lambda vertex:complex(np.mean(_tr(_dag(vertex)@R)))
    return v(W),v(Jy),v(Jx)

def spectral(s=Spec()):
    """Independent eigenstate band-pair sum on periodic k-grid."""
    D=_matrix_data(s)
    density,long,trans=_response(D['R'],D['W'],D['Jx'],D['Jy'])
    w=D['s'].omega+1j*D['s'].eta
    rhs=complex(np.mean(_tr(D['Fm'])-_tr(D['Fp'])))
    return {'density':density,'longitudinal':long,'transverse':trans,
            'q':D['q'],'ward_residual':complex(w*density-D['q']*long-rhs),
            'global_density_fsum':rhs,
            'electric_normalized_transverse':complex(trans/(-1j*D['q'])),
            'status':'FINITE_Q_PERIODIC_SPIN_ZEEMAN_SCALAR_DRIVE_NOT_SOURCE_PHE_DC'}

def liouville_independent(s=Spec()):
    """Independent dense 4x4 orbital Sylvester solve, no band-pair denominator."""
    hp,hm,zp,zm,W,Jx,Jy,q=_entries(s)
    ep,up=np.linalg.eigh(hp);em,um=np.linalg.eigh(hm)
    fp=fermi(ep,s.mu,s.temperature);fm=fermi(em,s.mu,s.temperature)
    Fp=(up*fp[...,None,:])@_dag(up)
    Fm=(um*fm[...,None,:])@_dag(um)
    src=W@Fm-Fp@W
    R=np.empty_like(src)
    w=s.omega+1j*s.eta
    for i in range(s.n):
        for j in range(s.n):
            A=np.kron(w*I-hp[i,j],I)+np.kron(I,hm[i,j].T)
            R[i,j]=np.linalg.solve(A,src[i,j].reshape(4)).reshape(2,2)
    den,long,trans=_response(R,W,Jx,Jy)
    return {'density':den,'longitudinal':long,'transverse':trans,
            'electric_normalized_transverse':trans/(-1j*q)}

def analytic_zeeman_mixed(s=Spec()):
    """Analytic B derivative, spin Z(k)=sigma_x only, NOT orbital-B current."""
    D=_matrix_data(s)
    ep,up=D['ep'],D['up'];em,um=D['em'],D['um']
    dFp=_frechet_f(ep,up,D['zp'],s)
    dFm=_frechet_f(em,um,D['zm'],s)
    dS=D['W']@dFm-dFp@D['W']
    Q=dS+D['zp']@D['R']-D['R']@D['zm']
    rb=(_dag(up)@Q@um)/D['den']
    dR=up@rb@_dag(um)
    den,lon,tra=_response(dR,D['W'],D['Jx'],D['Jy'])
    return {'dB_density':den,'dB_longitudinal':lon,'dB_transverse':tra,
            'dB_E_transverse':tra/(-1j*D['q']),
            'dB_ward_residual':(s.omega+1j*s.eta)*den-D['q']*lon,
            'status':'ANALYTIC_SPIN_B_DERIVATIVE_FINITE_Q_NO_ORBITAL_FIELD'}

def field_difference(s=Spec(),step=1e-4,points=5):
    if not np.isfinite(step) or step<=0 or points not in (3,5):
        raise ValueError('positive step and 3/5 point difference required')
    def f(B):return spectral(replace(s,field=B))
    x=s.field
    if points==3:
        a=f(x+step);b=f(x-step)
        return {k:(a[k]-b[k])/(2*step) for k in ('density','longitudinal','transverse','electric_normalized_transverse')}
    a=f(x-2*step);b=f(x-step);c=f(x+step);d=f(x+2*step)
    return {k:(a[k]-8*b[k]+8*c[k]-d[k])/(12*step) for k in ('density','longitudinal','transverse','electric_normalized_transverse')}

def incorrect_identity_vertex(s=Spec()):
    """Negative control: incorrectly replace U(k+q)^dagger U(k) by identity."""
    if abs(s.orbital_gauge)<1e-8:raise ValueError('requires momentum-dependent orbital basis')
    D=_matrix_data(s)
    wrong=np.broadcast_to(I,D['W'].shape)
    src=wrong@D['Fm']-D['Fp']@wrong
    rb=(_dag(D['up'])@src@D['um'])/D['den']
    R=D['up']@rb@_dag(D['um'])
    den,lon,trans=_response(R,wrong,D['Jx'],D['Jy'])
    return {'density':den,'longitudinal':lon,'transverse':trans,
            'ward_defect':(s.omega+1j*s.eta)*den-D['q']*lon}

def incorrect_local_velocity_vertex(s=Spec()):
    """Negative control: endpoint vy average violates exact finite-q Ward."""
    D=_matrix_data(s)
    k=2*np.pi*np.arange(s.n)/s.n
    _,y=np.meshgrid(k,k,indexing='ij')
    vy=(s.alpha*np.sin(y))[...,None,None]*I-np.cos(y)[...,None,None]*SX+(s.wilson*np.sin(y))[...,None,None]*SZ
    vymid=(vy+np.roll(vy,-s.q_steps,axis=1))*.5
    _,_,_,U=_field_and_vertices(s)
    Up=np.roll(U,-s.q_steps,axis=1)
    Jwrong=_dag(Up)@vymid@U
    wrong_long=complex(np.mean(_tr(_dag(Jwrong)@D['R'])))
    correct_density=spectral(s)['density']
    return {'longitudinal':wrong_long,
            'ward_defect':(s.omega+1j*s.eta)*correct_density-D['q']*wrong_long}
