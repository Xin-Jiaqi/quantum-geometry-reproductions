"""PT-degenerate staggered spin-geometry candidate: projector, U(2), Fermi-surface gates.
EXPLORATORY synthetic model; not the Xiang 2025 PRL reproduction or a new effect claim.
"""
import numpy as np
I=np.eye(2,dtype=complex); sx=np.array([[0,1],[1,0]],complex)
sy=np.array([[0,-1j],[1j,0]],complex); sz=np.diag([1,-1]).astype(complex)
K=np.kron; I4=np.eye(4,dtype=complex); L=K(sz,I); S=K(I,sy); A=K(sx,1j*sy)

def bloch(x,y,gamma=.31,b=.72,t=.65,alpha=.55,m=.72):
    e=-2*b*(np.cos(x)+np.cos(y))+gamma*np.sin(x)
    v0=2*b*np.sin(x)+gamma*np.cos(x)
    D=t*K(sx,I)+alpha*(np.sin(y)*K(sz,sx)-np.sin(x)*K(sz,sy))+m*K(sz,sz)
    r=np.sqrt(t*t+alpha**2*(np.sin(x)**2+np.sin(y)**2)+m*m)
    H=e*I4+D; v=v0*I4-alpha*np.cos(x)*K(sz,sy)
    Pm=(I4-D/r)/2; Pp=(I4+D/r)/2
    j=(v@S+S@v)/4; jm=(L@j+j@L)/2
    return H,v,Pm,Pp,j,jm,2*r

def sample(N=51,gamma=.31,b=.72,T=.15,mu=0):
    geom=np.zeros(2); surface=np.zeros(2); vmax=-1e10; cmin=1e10
    for x in np.linspace(-np.pi,np.pi,N,endpoint=False):
      for y in np.linspace(-np.pi,np.pi,N,endpoint=False):
        H,v,pm,pp,j,jm,d=bloch(x,y,gamma,b)
        geom+=np.array([np.trace(pm@j@pp@v).real, np.trace(pm@jm@pp@v).real])/d**2/N**2
        for P in (pm,pp):
          energy=np.trace(P@H).real/2; vel=np.trace(P@v).real/2
          if P is pm: vmax=max(vmax,energy)
          else: cmin=min(cmin,energy)
          minus_fprime=1/(4*T*np.cosh(np.clip((energy-mu)/(2*T),-60,60))**2)
          surface+=minus_fprime*vel*np.array([np.trace(P@j).real,np.trace(P@jm).real])/N**2
    return geom,surface,cmin-vmax

def checks():
    H,v,pm,pp,j,jm,d=bloch(.8,.6)
    assert np.linalg.norm(A@A.conj()+I4)<1e-12
    assert np.linalg.norm(A@H.conj()@A.conj().T-H)<1e-12
    assert np.linalg.norm(A@j.conj()@A.conj().T+j)<1e-12
    assert np.linalg.norm(A@jm.conj()@A.conj().T-jm)<1e-12
    assert np.linalg.norm(pm@pm-pm)<1e-12
    q=np.trace(pm@jm@pp@v).real/d**2
    E,U=np.linalg.eigh(H); rng=np.random.default_rng(241)
    def randU():
      Z=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2));Q,R=np.linalg.qr(Z)
      return Q@np.diag(np.diag(R)/abs(np.diag(R)))
    W=np.zeros((4,4),complex);W[:2,:2]=randU();W[2:,2:]=randU()
    Ug=U@W
    rotated=np.trace((Ug.conj().T@jm@Ug)[:2,2:]@(Ug.conj().T@v@Ug)[2:,:2]).real/d**2
    assert abs(q-rotated)<1e-12
    assert abs(np.trace(pm@j@pp@v).real/d**2)<1e-12
    return q

if __name__=='__main__':
    print('PT and U(2) checks; local staggered kernel',checks())
    for gamma in (0.,.15,.31,-.31):
      g,s,gap=sample(gamma=gamma)
      print('metal gamma',gamma,'integrated geometry',g,'Fermi-surface coefficient',s,'indirect gap',gap)
      assert abs(g[0])<1e-12 and abs(s[0])<1e-12
    g,s,gap=sample(gamma=.31,b=.12,T=.02)
    print('insulator gamma .31',g,s,'indirect gap',gap)
    assert gap>0 and abs(s[1])<1e-8
