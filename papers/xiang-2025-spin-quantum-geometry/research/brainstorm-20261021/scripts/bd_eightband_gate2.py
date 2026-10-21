"""B+D Gate 2: eight-band PT bilayer AC spin-density experiment.

Exploratory model only; not original Xiang/Jin/Wang PRL code, FE, or real DFT.
Two opposite-parity orbitals, two layers and spin. Natural units.
Requires numpy; run: python bd_eightband_gate2.py
"""
import json
import numpy as np

I=np.eye(2,dtype=complex)
X=np.array([[0,1],[1,0]],dtype=complex)
Y=np.array([[0,-1j],[1j,0]],dtype=complex)
Z=np.diag([1,-1]).astype(complex)
def kron(a,b,c): return np.kron(np.kron(a,b),c)
I8=np.eye(8,dtype=complex)
Tz=kron(Z,I,I); sx=kron(I,I,X); sy=kron(I,I,Y); sz=kron(I,I,Z)
Tx=kron(X,I,I); Ty=kron(Y,I,I)
PT=kron(X,Z,1j*Y); C2=kron(I,I,1j*Z)
Sminus_z=Tz@sz/2

def bond_sum(kx,ky,dx,dy,amplitude,axis=None,ell=.72,rmax=4):
    value=deriv=0j
    for nx in range(-rmax,rmax+1):
      for ny in range(-rmax,rmax+1):
        rx,ry=nx+dx,ny+dy
        weight=amplitude*np.exp(-(rx*rx+ry*ry)/(2*ell*ell))
        if axis=='x': weight*=rx
        elif axis=='y': weight*=ry
        piece=weight*np.exp(1j*(kx*rx+ky*ry))
        value+=piece;deriv+=1j*rx*piece
    return value,deriv

def interlayer(kx,ky,dx,dy,t0=.3,lam=.22):
    e,_=bond_sum(kx,ky,dx,dy,t0)
    x,_=bond_sum(kx,ky,dx,dy,lam,'x')
    y,_=bond_sum(kx,ky,dx,dy,lam,'y')
    off=lambda f,B:f.real*Tx@B-f.imag*Ty@B
    return (off(e,kron(I,I,I)+.18*kron(I,Z,I))
            +off(x,kron(I,X,Y))+off(y,kron(I,X,X)))

def model(kx,ky,dx,dy,t0=.3,lam=.22,orbital_gap=1.5,m=.48):
    e,de=bond_sum(kx,ky,dx,dy,t0)
    x,ex=bond_sum(kx,ky,dx,dy,lam,'x')
    y,ey=bond_sum(kx,ky,dx,dy,lam,'y')
    even=kron(I,I,I)+.18*kron(I,Z,I)
    oddx=kron(I,X,Y);oddy=kron(I,X,X)
    off=lambda f,B:f.real*Tx@B-f.imag*Ty@B
    H=(-.12*(np.cos(kx)+np.cos(ky))*I8
       +(orbital_gap+.17*(np.cos(kx)-np.cos(ky)))*kron(I,Z,I)
       +.30*Tz@(np.sin(ky)*sx-np.sin(kx)*sy)+m*Tz@sz
       +off(e,even)+off(x,oddx)+off(y,oddy))
    V=(.12*np.sin(kx)*I8-.17*np.sin(kx)*kron(I,Z,I)
       -.30*np.cos(kx)*Tz@sy
       +off(de,even)+off(ex,oddx)+off(ey,oddy))
    return H,V

def response(dx=.17,dy=.23,t0=.30,lam=.22,omega=.30,eta=.02,N=17,
             projector_check=False):
    """Complete occupied-to-empty *interband* insulator length-gauge Kubo.

    r_cv=-i v_cv/(E_c-E_v); chi=Q/(w-D+ieta)-Q*/(w+D+i eta).
    Charge and conventional overall signs are not fixed; only relative tests.
    """
    chi=0j; geom=0.; direct=np.inf;maxv=-np.inf;minc=np.inf
    projector_error=0.
    for kx in np.linspace(-np.pi,np.pi,N,endpoint=False):
      for ky in np.linspace(-np.pi,np.pi,N,endpoint=False):
        H,V=model(kx,ky,dx,dy,t0,lam)
        E,U=np.linalg.eigh(H)
        direct=min(direct,E[4]-E[3])
        maxv=max(maxv,E[3]);minc=min(minc,E[4])
        S=U.conj().T@Sminus_z@U
        v=U.conj().T@V@U
        local=0j
        for n in range(4):
          for m in range(4,8):
            gap=E[m]-E[n]
            q=S[n,m]*v[m,n]
            a=-1j*q/gap
            local+=a/(omega-gap+1j*eta)-a.conjugate()/(omega+gap+1j*eta)
            geom+=q.real/gap**2/N**2
        chi+=local/N**2
        if projector_check:
          proj=0j
          for n in (0,2):
            p=U[:,n:n+2]@U[:,n:n+2].conj().T
            for m in (4,6):
              q=U[:,m:m+2]@U[:,m:m+2].conj().T
              d=E[m:m+2].mean()-E[n:n+2].mean()
              a=-1j*np.trace(p@Sminus_z@q@V)/d
              proj+=a/(omega-d+1j*eta)-a.conjugate()/(omega+d+1j*eta)
          projector_error=max(projector_error,abs(proj-local))
    return chi,geom,direct,minc-maxv,projector_error

def tests():
    assert np.linalg.norm(PT@PT.conj()+I8)<1e-12
    for dx,dy in [(0,0),(.17,.23),(-.17,-.23),(.32,-.11)]:
      for kx,ky in [(0.3,-.8),(.84,.53),(-1.51,.68)]:
        H,V=model(kx,ky,dx,dy)
        Hp,Vp=model(-kx,-ky,-dx,-dy)
        assert np.linalg.norm(H-H.conj().T)<5e-13
        assert np.linalg.norm(PT@H.conj()@PT.conj().T-H)<5e-13
        assert np.linalg.norm(PT@V.conj()@PT.conj().T-V)<5e-13
        assert np.linalg.norm(C2@H@C2.conj().T-Hp)<5e-13
        assert np.linalg.norm(C2@V@C2.conj().T+Vp)<5e-13
        W=kron(np.diag(np.exp(1j*np.pi*np.array([dx,-dx]))),I,I)
        Hx,Vx=model(kx+2*np.pi,ky,dx,dy)
        assert np.linalg.norm(W@H@W.conj().T-Hx)<5e-13
        assert np.linalg.norm(W@V@W.conj().T-Vx)<5e-13
        hi=interlayer(kx,ky,dx,dy)
        torque=lambda a:1j*(a@Sminus_z-Sminus_z@a)
        assert np.linalg.norm(torque(H)-torque(hi)-torque(H-hi))<5e-13
        e=np.linalg.eigvalsh(H)
        assert max(abs(e[0::2]-e[1::2]))<1e-12
    for B in (kron(X,X,X),kron(X,X,Y),kron(Y,X,X),kron(Y,X,Y)):
      assert np.linalg.norm(PT@B.conj()@PT.conj().T-B)<5e-13
    rng=np.random.default_rng(731)
    H,V=model(.62,.41,.17,.23)
    E,U=np.linalg.eigh(H)
    W=np.zeros((8,8),complex)
    for i in range(0,8,2):
      a=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
      q,r=np.linalg.qr(a);W[i:i+2,i:i+2]=q@np.diag(np.diag(r)/abs(np.diag(r)))
    def invariant(u):
      a=u.conj().T@Sminus_z@u;b=u.conj().T@V@u
      return sum((a[n,m]*b[m,n]).real/(E[m]-E[n])**2
                 for n in range(4) for m in range(4,8))
    assert abs(invariant(U)-invariant(U@W))<1e-12
    a,_=bond_sum(.49,.51,.17,.23,.3,rmax=6)
    b,_=bond_sum(.49,.51,1.17,.23,.3,rmax=6)
    assert abs(a-b)<1e-9

def run():
    tests()
    data=[]
    for N in (11,17,25):
      for name,(t0,lam) in dict(none=(0,0),ordinary_only=(.30,0),
                                spinful_only=(0,.22),both=(.30,.22)).items():
        for sign in (1,-1) if name=='both' else (1,):
          c,g,dr,ind,err=response(sign*.17,sign*.23,t0,lam,N=N,
                     projector_check=(N==11 and name=='both'))
          assert ind>.60 and dr>.80
          data.append(dict(N=N,case=name,sign=sign,chi_real=c.real,
                           chi_imag=c.imag,geometry_kernel=g,
                           direct_gap=dr,indirect_gap=ind,projector_err=err))
    for N in (11,17,25):
      rows=[d for d in data if d['N']==N]
      a=next(d for d in rows if d['case']=='both' and d['sign']==1)
      b=next(d for d in rows if d['case']=='both' and d['sign']==-1)
      c=complex(a['chi_real'],a['chi_imag'])
      d=complex(b['chi_real'],b['chi_imag'])
      assert abs(c+d)<5e-13
      assert abs(c.imag)>2e-7
      assert all(abs(complex(z['chi_real'],z['chi_imag']))<2e-12
                 for z in rows if z['case']!='both')
      print('N=%d: chi(+delta)=%.11g%+.11gj, indirect gap %.8f, projector error %.2g'%
            (N,c.real,c.imag,a['indirect_gap'],a['projector_err']))
    z=response(dx=0,dy=0,N=17)[0]
    assert abs(z)<1e-12
    print('PASS: symmetry, U2, projectors, switch parity, dual-channel and insulator gates')
    print(json.dumps(dict(parameters=dict(omega=.30,eta=.02,delta=[.17,.23],
                                         ordinary_t0=.30,spinful_lambda=.22),
                          data=data),indent=2))
if __name__=='__main__':run()
