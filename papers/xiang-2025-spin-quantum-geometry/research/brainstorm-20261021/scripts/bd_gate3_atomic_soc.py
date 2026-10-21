"""B+D Gate 3: atomic SOC -> interlayer spin hopping; competing-response audit.
Independent research model, NOT the Xiang PRL reproduction, NOT a real
ferroelectric material, and NOT a new-effect claim.

Run: python scripts/bd_gate3_atomic_soc.py (numpy only).
Orbitals (s,px,py,pz) per layer, spin; full size 16.
"""
import numpy as np

I=np.eye(2,dtype=complex)
X=np.array([[0,1],[1,0]],complex)
Y=np.array([[0,-1j],[1j,0]],complex)
Z=np.diag([1,-1]).astype(complex)
K=lambda a,b,c: np.kron(np.kron(a,b),c)
DP=np.diag([1,-1,-1,-1]).astype(complex)
DC2=np.diag([1,-1,-1,1]).astype(complex)
DMZ=np.diag([1,1,1,-1]).astype(complex)
PT=K(X,DP,1j*Y)
C2=K(I,DC2,1j*Z)
MZT=K(X,DMZ,X)
TR=K(I,np.eye(4),1j*Y)
TZ=K(Z,np.eye(4),I)
S_TOT_Z=K(I,np.eye(4),Z)/2
S_STAG_Z=TZ@S_TOT_Z
EORB=np.diag([-1.15,4.,4.,1.15]).astype(complex)
LX=np.zeros((4,4),complex); LY=LX.copy(); LZ=LX.copy()
LX[2,3]=-1j; LX[3,2]=1j
LY[1,3]=1j; LY[3,1]=-1j
LZ[1,2]=-1j; LZ[2,1]=1j
LOW=np.array([8*a+2*b+c for a in range(2) for b in (0,3) for c in range(2)])
HIGH=np.array([i for i in range(16) if i not in LOW])

def model(kx,ky,dx=.17,dy=.23,t_scale=1.,v_sp=.23,lam_soc=.36,
          vertical=.95,m=.36,ell=.68,rmax=3,polar_u=0.):
    """Atomic real-space, spin-INDEPENDENT interlayer Slater-Koster bonds.
    All spin-orbit coupling occurs on the same p-atom via L dot S.
    polar_u is an EXTERNAL potential proxy, not self-consistent FE.
    """
    onsite_orb=EORB+np.diag([-.12*(np.cos(kx)+np.cos(ky)),
                            0,0,.12*(np.cos(kx)+np.cos(ky))])
    H=K(I,onsite_orb,I)+m*K(Z,np.eye(4),Z)+polar_u*TZ
    H+=lam_soc/2*(K(I,LX,X)+K(I,LY,Y)+K(I,LZ,Z))
    V=K(I,np.diag([.12*np.sin(kx),0,0,-.12*np.sin(kx)]),I)
    t=np.zeros((4,4),complex);tv=t.copy()
    for nx in range(-rmax,rmax+1):
      for ny in range(-rmax,rmax+1):
        rx,ry=nx+dx,ny+dy
        distance=np.sqrt(rx*rx+ry*ry+vertical**2)
        direction=np.array([rx,ry,vertical])/distance
        B=np.zeros((4,4))
        B[0,0]=.28
        B[0,1:]=v_sp*direction
        B[1:,0]=-v_sp*direction
        B[1:,1:]=.04*np.eye(3)+.17*np.outer(direction,direction)
        a=t_scale*np.exp(-(rx*rx+ry*ry)/(2*ell*ell))
        piece=a*np.exp(1j*(kx*rx+ky*ry))
        t+=piece*B;tv+=1j*rx*piece*B
    T=np.kron(t,I);Tv=np.kron(tv,I)
    H[:8,8:]+=T;H[8:,:8]+=T.conj().T
    V[:8,8:]+=Tv;V[8:,:8]+=Tv.conj().T
    return H,V

def downfold(kx=.5,ky=.4,dx=.17,dy=.23,v_sp=.23,lam_soc=.36):
    """E=0 Schur complement, a LOW-energy approximation, not exact SW."""
    H,_=model(kx,ky,dx,dy,v_sp=v_sp,lam_soc=lam_soc)
    ll=H[np.ix_(LOW,LOW)];lh=H[np.ix_(LOW,HIGH)]
    hh=H[np.ix_(HIGH,HIGH)]
    heff=ll-lh@np.linalg.solve(hh,lh.conj().T)
    T=heff[:4,4:]
    s_to_pz=T[:2,2:4]
    coeff=np.array([np.trace(P@s_to_pz)/2 for P in (I,X,Y,Z)])
    return coeff,heff,ll

def susceptibility(N=13,dx=.17,dy=.23,m=.36,u=0.,t_scale=1.,
                   v_sp=.23,lam_soc=.36,omega=.3,eta=.02,
                   observable=S_STAG_Z):
    """All occupied (4) / empty (12) length-gauge interband Kubo pairs.
    Natural units, charge and overall response sign convention not SI.
    """
    ret=0j;maxv=-np.inf;minc=np.inf;mind=np.inf;layer=0.
    for kx in np.linspace(-np.pi,np.pi,N,endpoint=False):
      for ky in np.linspace(-np.pi,np.pi,N,endpoint=False):
        H,V=model(kx,ky,dx,dy,t_scale,v_sp,lam_soc,
                  m=m,polar_u=u)
        E,U=np.linalg.eigh(H);maxv=max(maxv,E[3]);minc=min(minc,E[4])
        mind=min(mind,E[4]-E[3])
        P=U[:,:4]@U[:,:4].conj().T
        layer+=np.trace(P@TZ).real/N**2
        A=U.conj().T@observable@U; W=U.conj().T@V@U
        d=E[4:][None,:]-E[:4,None]
        q=A[:4,4:]*W[4:,:4].T
        a=-1j*q/d
        ret+=np.sum(a/(omega-d+1j*eta)-a.conj()/(omega+d+1j*eta))/N**2
    return ret,mind,minc-maxv,layer

def occupied_proj(kx,ky,dx,dy,m,bz=0.):
    H,_=model(kx,ky,dx,dy,m=m)
    H=H+bz*S_TOT_Z
    _,U=np.linalg.eigh(H)
    return U[:,:4]@U[:,:4].conj().T

def mixed_curvature(N=11,m=.36,dx=.17,dy=.23,h=.001):
    """Independent finite-difference occupied projector on (kx, Bz)."""
    val=0.
    for kx in np.linspace(-np.pi,np.pi,N,endpoint=False):
      for ky in np.linspace(-np.pi,np.pi,N,endpoint=False):
        P=occupied_proj(kx,ky,dx,dy,m)
        dk=(occupied_proj(kx+h,ky,dx,dy,m)-
            occupied_proj(kx-h,ky,dx,dy,m))/(2*h)
        db=(occupied_proj(kx,ky,dx,dy,m,h)-
            occupied_proj(kx,ky,dx,dy,m,-h))/(2*h)
        val+=(1j*np.trace(P@(dk@db-db@dk))).real/N**2
    return val

def check_symmetries():
    assert np.linalg.norm(PT@PT.conj()+np.eye(16))<1e-12
    for kx,ky in ((.31,.45),(-.55,.93)):
      for dx,dy in ((0,0),(.17,.23),(-.17,-.23)):
        H,V=model(kx,ky,dx,dy)
        HP,VP=model(-kx,-ky,-dx,-dy)
        assert np.linalg.norm(H-H.conj().T)<1e-12
        assert np.linalg.norm(PT@H.conj()@PT.conj().T-H)<1e-12
        assert np.linalg.norm(PT@V.conj()@PT.conj().T-V)<1e-12
        assert np.linalg.norm(C2@H@C2.conj().T-HP)<1e-12
        assert np.linalg.norm(C2@V@C2.conj().T+VP)<1e-12
        hu,vu=model(kx,ky,dx,dy,polar_u=.08)
        hd,vd=model(-kx,-ky,-dx,-dy,polar_u=-.08)
        assert np.linalg.norm(MZT@hu.conj()@MZT.conj().T-hd)<1e-12
        assert np.linalg.norm(MZT@vu.conj()@MZT.conj().T+vd)<1e-12
        assert np.linalg.norm(PT@hu.conj()@PT.conj().T-hu)>.1
        assert np.linalg.norm(C2@hu@C2.conj().T-hd)>.1
        he,_,_=downfold(kx,ky,dx,dy)
        pt8=PT[np.ix_(LOW,LOW)]
        assert np.linalg.norm(pt8@he.conj()@pt8.conj().T-he)<1e-12
    H,_=model(.4,.5,m=0.)
    HR,_=model(-.4,-.5,m=0.)
    assert np.linalg.norm(TR@H.conj()@TR.conj().T-HR)<1e-12

def audit():
    check_symmetries()
    print("PASS: full 16-band PT, C2z, m_z*T daughter connection, TR at m=0")
    for vsp,soc in ((.23,.36),(.23,0.),(0.,.36)):
      c,_,_=downfold(v_sp=vsp,lam_soc=soc)
      print("downfold s->pz (I,sigma_x,sigma_y,sigma_z):",
            "v_sp",vsp,"atomic SOC",soc,np.array2string(c,precision=8))
      if vsp==0 or soc==0: assert np.linalg.norm(c[1:])<1e-11
      else: assert np.linalg.norm(c[1:])>1e-3
    for m in (0.,.36,-.36):
      a,dr,ind,_=susceptibility(m=m,observable=S_STAG_Z)
      b,_,_,_=susceptibility(m=m,observable=S_TOT_Z)
      print("m",m,"indgap",ind,"stag_z(AC)",a,"total_z(AC)",b)
      assert ind>0.9 and dr>1.
    for kwargs in ({"t_scale":0.},{"v_sp":0.},{"lam_soc":0.}):
      a,_,_,_=susceptibility(observable=S_TOT_Z,**kwargs)
      assert abs(a)<1e-11
    chi0,_,_,_=susceptibility(N=11,omega=0.,eta=0.,observable=S_TOT_Z)
    mixed=mixed_curvature(N=11)
    print("magnetic static chi_total_z",chi0,
          "mixed F(kx,Bz)",mixed,"sum",chi0.real+mixed)
    assert abs(chi0.real+mixed)<3e-11
    for sign in (1,-1):
      dx,dy,u=sign*.17,sign*.23,sign*.08
      total,_,_,layer=susceptibility(N=13,dx=dx,dy=dy,u=u,
                          omega=0.,eta=0.,observable=S_TOT_Z)
      stag,_,_,_=susceptibility(N=13,dx=dx,dy=dy,u=u,
                          omega=0.,eta=0.,observable=S_STAG_Z)
      print("daughter",sign,"layer occupancy difference",layer,
            "static total_z",total,"static staggered_z",stag)
    print("PASS: micro-SOC, magnetic vs nonmagnetic, mixed-curvature",
          "and polar-proxy controls. Novelty: NOT ESTABLISHED.")

if __name__=="__main__":
    audit()
