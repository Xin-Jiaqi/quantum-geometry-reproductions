"""B+D Gate 4A: explicit two-atom honeycomb stacking, toy energy and layer dipole.

One spinless p_z-like orbital at each inequivalent A/B site on each layer.
A/B chemistry is model-only (hBN-inspired), NOT a magnetic hBN material.
The total toy energy is an independent band-energy plus screened ionic pair
interaction. NOT self-consistent DFT, SI polarization, or an ab initio barrier.
Run directly: python scripts/bd_gate4_polar_honeycomb.py (numpy only).
"""
import numpy as np

A1=np.array([1.,0.]);A2=np.array([.5,np.sqrt(3)/2])
B=(A1+A2)/3
AREA=np.linalg.det(np.column_stack((A1,A2)))
B1=np.array([2*np.pi,-2*np.pi/np.sqrt(3)])
B2=np.array([0.,4*np.pi/np.sqrt(3)])
ORBS=(np.zeros(2),B)

def fraction_delta(frac):return frac[0]*A1+frac[1]*A2

def rvecs(rmax):
    return np.array([a*A1+b*A2
                     for a in range(-rmax,rmax+1)
                     for b in range(-rmax,rmax+1)])

def model(k,delta,delta_chem=1.4,t_in=.65,t_out=.35,ell=.33,rmax=3):
    intra=-t_in*sum(np.exp(1j*k@d) for d in (B,B-A1,B-A2))
    h=np.array([[delta_chem,intra],[intra.conjugate(),-delta_chem]])
    T=np.zeros((2,2),complex)
    R=rvecs(rmax)
    for a in range(2):
        for b in range(2):
            d=R+delta+ORBS[b]-ORBS[a]
            amp=t_out*np.exp(-np.sum(d*d,axis=1)/(2*ell**2))
            T[a,b]=np.dot(amp,np.exp(1j*d@k))
    return np.block([[h,T],[T.conj().T,h]])

def bands(frac,N=11,delta_chem=1.4,t_out=.35):
    delta=fraction_delta(frac)
    imbalance=0.;total=0.;direct=1e9;vmax=-1e9;cmin=1e9
    for i in range(N):
        for j in range(N):
            k=(i*B1+j*B2)/N
            e,U=np.linalg.eigh(model(k,delta,delta_chem=delta_chem,t_out=t_out))
            P=U[:,:2]@U[:,:2].conj().T
            imbalance+=(np.trace(P[:2,:2]).real-np.trace(P[2:,2:]).real)/N**2
            total+=e[:2].sum()/N**2
            direct=min(direct,e[2]-e[1])
            vmax=max(vmax,e[1]);cmin=min(cmin,e[2])
    return np.array([imbalance,total,direct,cmin-vmax])

def ion_pair(frac,charge=.9,screen=2.,core=.12,core_w=.25,h=.75,rmax=7):
    delta=fraction_delta(frac);R=rvecs(rmax);value=0.
    for a in range(2):
        for b in range(2):
            d=R+delta+ORBS[b]-ORBS[a]
            r=np.sqrt(np.sum(d*d,axis=1)+h*h)
            value+=np.sum(charge**2*(1 if a==b else -1)*
               np.exp(-screen*r)/r+core*np.exp(-(r-h)**2/(2*core_w**2)))
    return value

def energy(frac,N=11,charge=.9):
    return bands(frac,N)[1]+ion_pair(frac,charge=charge)

def hessian(frac,h=.001,N=13,charge=.9):
    f=np.array(frac);out=np.zeros((2,2))
    for a in range(2):
        for b in range(2):
            ea=np.eye(2)[a]*h;eb=np.eye(2)[b]*h
            out[a,b]=(energy(f+ea+eb,N,charge)-energy(f+ea-eb,N,charge)
                     -energy(f-ea+eb,N,charge)+energy(f-ea-eb,N,charge))/(4*h*h)
    return out

def audit():
    f=(1/3,1/3);opp=(-1/3,-1/3);cent=(0.,0.)
    for N in (7,11,19):
        a,b,c=bands(f,N),bands(opp,N),bands(cent,N)
        assert abs(a[0]+b[0])<1e-12 and abs(a[1]-b[1])<1e-12
        assert abs(c[0])<1e-12 and a[3]>2.7
        assert energy(cent,N)>energy(f,N)
        print("N",N,"layer electron imbalance ±AB/BA",a[0],b[0],
              "toy E_AB",energy(f,N),"E_AA",energy(cent,N),"indirect gap",a[3])
    assert abs(bands(f,19,t_out=0)[0])<1e-12
    assert abs(bands(f,19,delta_chem=0)[0])<1e-12
    minima=hessian(f)
    opposite=hessian(opp)
    saddle=hessian((0.,.5))
    assert np.min(np.linalg.eigvalsh(minima))>.10
    assert np.min(np.linalg.eigvalsh(opposite))>.10
    ev=np.linalg.eigvalsh(saddle)
    assert ev[0]<0<ev[1]
    print("AB Hessian",np.linalg.eigvalsh(minima))
    print("candidate saddle (0,1/2) Hessian",ev,
          "saddle-minus-AB toy energy",energy((0.,.5),19)-energy(f,19))
    print("AA-minus-AB toy energy",energy(cent,19)-energy(f,19))
    # Baseline local minima depend on ionic parameters: q=0,.45 do not pass.
    for q in (0.,.45,.9,1.3):
        print("ionic q",q,"AB Hessian",np.linalg.eigvalsh(hessian(f,charge=q)))
    exchange=np.kron(np.array([[0.,1.],[1.,0.]]),np.eye(2))
    for k in (np.array([.1,.2]),np.array([1.2,-.4])):
        for ff in ((.19,.24),f):
            d=fraction_delta(ff)
            assert np.linalg.norm(exchange@model(k,d)@exchange-model(k,-d))<1e-12
    for ff in ((.17,.11),f):
        base=ion_pair(ff,rmax=7)
        for lattice in ((1,0),(0,1)):
            assert abs(ion_pair(np.array(ff)+lattice,rmax=7)-base)<2e-5
    print("PASS: mirror pairing, band gap, polar charge, Hessian, disconnected and identical-site controls")

if __name__=="__main__":audit()
