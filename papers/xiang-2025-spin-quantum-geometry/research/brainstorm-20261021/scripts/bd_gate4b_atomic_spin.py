"""Gate 4B: 32-band honeycomb bilayer with atom-local SOC and fixed A-type AFM.

Physical A/B basis geometry, s, px, py, pz orbitals, two layers, spin.
All interlayer hoppings are spin-INDEPENDENT real-space Slater-Koster bonds.
Spin mixing derives exclusively from on-atom p-orbital L dot S.
This is a SYNTHETIC magnetic honeycomb, NOT magnetic hBN, DFT or a
novel response. Finite-field static magnetoelectric Kubo, not driven DC.

Requires bd_gate4_polar_honeycomb.py in the same directory and numpy.
Run: python scripts/bd_gate4b_atomic_spin.py
"""
import numpy as np
from bd_gate4_polar_honeycomb import (A1,A2,B,B1,B2,ORBS,rvecs,
                                      fraction_delta,ion_pair)
I=np.eye(2,dtype=complex)
X=np.array([[0,1],[1,0]],complex)
Y=np.array([[0,-1j],[1j,0]],complex)
Z=np.diag([1,-1]).astype(complex)
I4=np.eye(4,dtype=complex)
K=lambda a,b,c,d:np.kron(np.kron(np.kron(a,b),c),d)
LX=np.zeros((4,4),complex);LY=LX.copy();LZ=LX.copy()
LX[2,3]=-1j;LX[3,2]=1j
LY[1,3]=1j;LY[3,1]=-1j
LZ[1,2]=-1j;LZ[2,1]=1j
Tspin_z=K(I,I,I4,Z)/2
STAG_z=K(Z,I,I4,Z)/2
Tspin_y=K(I,I,I4,Y)/2
STAG_y=K(Z,I,I4,Y)/2
Tspin_x=K(I,I,I4,X)/2
STAG_x=K(Z,I,I4,X)/2
MZT=K(X,I,np.diag([1,1,1,-1]),X)
TR=K(I,I,I4,1j*Y)
LAYER=K(Z,I,I4,I)

def sk(d,ss=.35,sp=.20,ppi=.03,psig=.16):
    n=np.asarray(d)/np.linalg.norm(d)
    M=np.zeros((4,4))
    M[0,0]=ss
    M[0,1:]=sp*n;M[1:,0]=-sp*n
    M[1:,1:]=ppi*np.eye(3)+psig*np.outer(n,n)
    return M

def model(k,delta,m=.34,soc=.30,out_scale=1.,h=.75,ell=.33,rc=3):
    h16=np.zeros((16,16),complex)
    v16=np.zeros((16,16),complex)
    def idx(layer,site):
        return slice(8*layer+4*site,8*layer+4*site+4)
    for l in range(2):
        for s in range(2):
            e=-.35 if s==0 else -2.75
            h16[idx(l,s),idx(l,s)]+=np.diag([e,4.1,4.1,1.65])
        for d in (B,B-A1,B-A2):
            term=sk([d[0],d[1],0],ss=-.62,sp=.24,
                    ppi=.06,psig=.11)*np.exp(1j*d@k)
            h16[idx(l,0),idx(l,1)]+=term
            h16[idx(l,1),idx(l,0)]+=term.conj().T
            v16[idx(l,0),idx(l,1)]+=1j*d[0]*term
            v16[idx(l,1),idx(l,0)]+=(1j*d[0]*term).conj().T
    R=rvecs(rc)
    for a in range(2):
        for b in range(2):
            bonds=R+delta+ORBS[b]-ORBS[a]
            for d in bonds:
                term=(out_scale*np.exp(-d@d/(2*ell*ell))*
                    sk([d[0],d[1],h])*np.exp(1j*d@k))
                h16[idx(0,a),idx(1,b)]+=term
                h16[idx(1,b),idx(0,a)]+=term.conj().T
                v16[idx(0,a),idx(1,b)]+=1j*d[0]*term
                v16[idx(1,b),idx(0,a)]+=(1j*d[0]*term).conj().T
    H=np.kron(h16,I);V=np.kron(v16,I)
    H+=m*K(Z,I,I4,Z)
    H+=soc/2*(K(I,I,LX,X)+K(I,I,LY,Y)+K(I,I,LZ,Z))
    return H,V

def response(frac,N=9,m=.34,soc=.30,out_scale=1.,omega=0.,eta=0.,
             S=Tspin_z,occ_count=8):
    delta=fraction_delta(frac);chi=0j;direct=1e10
    dipole_electrons=0.;vmax=-1e10;cmin=1e10
    for i in range(N):
        for j in range(N):
            k=(i*B1+j*B2)/N
            H,V=model(k,delta,m,soc,out_scale)
            E,U=np.linalg.eigh(H)
            P=U[:,:occ_count]@U[:,:occ_count].conj().T
            dipole_electrons+=np.trace(P@LAYER).real/N**2
            direct=min(direct,E[occ_count]-E[occ_count-1])
            vmax=max(vmax,E[occ_count-1])
            cmin=min(cmin,E[occ_count])
            O=U.conj().T@S@U
            W=U.conj().T@V@U
            d=E[occ_count:][None,:]-E[:occ_count,None]
            q=O[:occ_count,occ_count:]*W[occ_count:,:occ_count].T
            r=-1j*q/d
            chi+=np.sum(r/(omega-d+1j*eta)-
                        r.conj()/(omega+d+1j*eta))/N**2
    return chi,dipole_electrons,direct,cmin-vmax

def toy_total_energy(frac,N=5):
    """NOT self-consistent DFT. Band-energy plus independent pair potential."""
    delta=fraction_delta(frac);v=0.
    for i in range(N):
        for j in range(N):
            k=(i*B1+j*B2)/N
            E=np.linalg.eigvalsh(model(k,delta)[0])
            v+=E[:8].sum()/N**2
    return v+ion_pair(frac)

def audit(dense=False):
    for frac in ((1/3,1/3),(-1/3,-1/3),(.25,.13),(0,0)):
        d=fraction_delta(frac)
        for k in (np.array([.31,.83]),np.array([-1.2,.73])):
            H,V=model(k,d)
            Hr,Vr=model(-k,-d)
            assert np.linalg.norm(H-H.conj().T)<1e-12
            assert np.linalg.norm(MZT@H.conj()@MZT.conj().T-Hr)<1e-12
            assert np.linalg.norm(MZT@V.conj()@MZT.conj().T+Vr)<1e-12
    print("PASS: hermiticity and antiunitary m_z*T connection")
    for N in ((7,11,15,19,23) if dense else (7,11)):
        rows=[]
        for obs,label in ((Tspin_x,'total_x'),(STAG_x,'staggered_x'),
                          (Tspin_z,'total_z'),(STAG_z,'staggered_z')):
            pos=response((1/3,1/3),N=N,S=obs)
            neg=response((-1/3,-1/3),N=N,S=obs)
            par=+1 if label.startswith('total') else -1
            assert abs(pos[0]-par*neg[0])<1e-11
            assert pos[3]>0.1 and neg[3]>0.1
            if label.endswith('z'):
                assert abs(pos[0])<1e-12
            rows.append((label,pos[0].real,neg[0].real))
        print("mesh",N,"responses",rows)
    N=7
    zero_m=response((1/3,1/3),N=N,m=0,S=Tspin_x)[0]
    zero_soc=response((1/3,1/3),N=N,soc=0,S=Tspin_x)[0]
    zero_t=response((1/3,1/3),N=N,out_scale=0,S=Tspin_x)[0]
    assert max(abs(zero_m),abs(zero_soc),abs(zero_t))<1e-12
    chi,polar,_,gap=response((1/3,1/3),N=N,S=Tspin_x)
    assert abs(chi)>1e-4 and abs(polar)>.03 and gap>.10
    print("PASS: spin response disappears for m=0, SOC=0, or t_inter=0")
    print("32-band toy energy: AB",toy_total_energy((1/3,1/3)),
          "AA",toy_total_energy((0,0)),
          "(AB is a local, NOT global, minimum; check Hessian separately)")
if __name__=="__main__":
    import sys
    audit(dense='--dense' in sys.argv)
