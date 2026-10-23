"""Round26 in-memory same-Wannier-gauge physical-operator input contract.

Fourier X(q)=sum_R exp(2pi i R.q) X(R), where q is fractional reciprocal.
H [eV], projected physical Pauli spin [dimensionless], projected position
connection A [angstrom]. A includes orbital centers exactly once.
D_k X = d_k X-i[A,X]. No Wannier file importer/material/KME claim.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

UNITS={"energy":"eV","length":"angstrom","spin":"dimensionless_pauli"}

@dataclass
class WannierOperatorBundle:
    translations: np.ndarray
    lattice: np.ndarray
    h_r: np.ndarray
    spin_r: np.ndarray
    position_r: np.ndarray
    centers: np.ndarray
    provenance: dict
    units: dict

    def validate(self,tol=1e-9):
        t=np.asarray(self.translations)
        a=np.asarray(self.lattice,float)
        h=np.asarray(self.h_r,complex)
        s=np.asarray(self.spin_r,complex)
        p=np.asarray(self.position_r,complex)
        c=np.asarray(self.centers,float)
        if self.units!=UNITS:
            raise ValueError("explicit eV, angstrom and dimensionless_pauli units required")
        if (t.ndim!=2 or t.shape[1]!=2 or len(t)==0 or
            not np.all(np.isfinite(t)) or np.any(t!=np.rint(t))):
            raise ValueError("translations must be integer (NR,2)")
        t=t.astype(int)
        if len(set(map(tuple,t)))!=len(t):
            raise ValueError("duplicate translation")
        if a.shape!=(2,3) or not np.all(np.isfinite(a)) or np.linalg.matrix_rank(a)<2:
            raise ValueError("two independent lattice vectors required")
        if h.ndim!=3 or h.shape[0]!=len(t) or h.shape[1]!=h.shape[2] or h.shape[-1]<2:
            raise ValueError("h_r must be (NR,N,N), N>=2")
        n=h.shape[-1]
        if s.shape!=(len(t),3,n,n) or p.shape!=(len(t),3,n,n) or c.shape!=(n,3):
            raise ValueError("spin_r, position_r and centers shapes mismatch")
        if not all(np.all(np.isfinite(x)) for x in (h,s,p,c)):
            raise ValueError("nonfinite operator")
        meta=self.provenance
        if not isinstance(meta,dict):raise ValueError("operator provenance missing")
        ids=meta.get("basis_ids")
        if (not isinstance(ids,dict) or set(ids)!={"hamiltonian","spin","position"} or
            any(not isinstance(v,str) or not v.strip() for v in ids.values()) or
            len(set(ids.values()))!=1):
            raise ValueError("same Wannier basis ID required for H/spin/position")
        if meta.get("spin_kind")!="projected_physical_pauli":
            raise ValueError("physical Pauli spin required, not orbital pseudospin")
        if not isinstance(meta.get("spin_evidence"),str) or not meta["spin_evidence"].strip():
            raise ValueError("explicit physical spin evidence required")
        if meta.get("orthonormal_basis") is not True:
            raise ValueError("orthonormal Wannier basis declaration required")
        if meta.get("position_convention")!="full_connection_including_centers":
            raise ValueError("position connection must include centers exactly once")
        ids_R={tuple(r):i for i,r in enumerate(t)}
        if (0,0) not in ids_R:raise ValueError("missing onsite R=0")
        for i,r in enumerate(t):
            j=ids_R.get(tuple(-r))
            if j is None:raise ValueError("missing inverse -R")
            for x in (h,s,p):
                if not np.allclose(x[i],np.swapaxes(x[j].conj(),-1,-2),atol=tol,rtol=0):
                    raise ValueError("R/-R Hermiticity violation")
        actual=np.diagonal(p[ids_R[(0,0)]],axis1=-2,axis2=-1).T.real
        if not np.allclose(actual,c,atol=tol,rtol=0):
            raise ValueError("centers must be diag of onsite position connection")
        for op in self.fourier([0,0],validate=False)["spin"]:
            if np.linalg.norm(op,ord=2)>1+tol:
                raise ValueError("projected Pauli norm exceeds one")
        return self

    def fourier(self,q,*,validate=True):
        if validate:self.validate()
        q=np.asarray(q,float)
        if q.shape!=(2,) or not np.all(np.isfinite(q)):
            raise ValueError("finite fractional k shape (2,) required")
        t=np.asarray(self.translations,int)
        R_cart=t@np.asarray(self.lattice,float)
        phase=np.exp(2j*np.pi*(t@q))
        f=lambda x:np.einsum("r,r...->...",phase,x,optimize=True)
        H=f(np.asarray(self.h_r,complex))
        S=f(np.asarray(self.spin_r,complex))
        A=f(np.asarray(self.position_r,complex))
        dH=np.einsum("r,ra,rij->aij",phase,1j*R_cart,np.asarray(self.h_r,complex),optimize=True)
        dS=np.einsum("r,ra,rbij->abij",phase,1j*R_cart,np.asarray(self.spin_r,complex),optimize=True)
        DH=dH-1j*(A@H-H@A)
        DS=dS-1j*(A[:,None]@S[None]-S[None]@A[:,None])
        return dict(H=H,spin=S,connection=A,partial_H=dH,partial_spin=dS,
                    covariant_H=DH,covariant_spin=DS,scope="model/interface only; no material KME")

def global_gauge_transform(b,U,*,suffix="U"):
    b.validate();U=np.asarray(U,complex);n=b.h_r.shape[-1]
    if (U.shape!=(n,n) or not np.all(np.isfinite(U)) or
        not np.allclose(U.conj().T@U,np.eye(n),atol=1e-10)):
        raise ValueError("unitary (N,N) required")
    rot=lambda x:U.conj().T@x@U
    h,s,p=rot(b.h_r),rot(b.spin_r),rot(b.position_r)
    iz=np.flatnonzero(np.all(np.asarray(b.translations)==0,axis=1))[0]
    centers=np.diagonal(p[iz],axis1=-2,axis2=-1).T.real
    meta=dict(b.provenance)
    old=next(iter(meta["basis_ids"].values()))
    meta["basis_ids"]={x:old+"/"+suffix for x in ("hamiltonian","spin","position")}
    meta["spin_evidence"]=meta["spin_evidence"]+"; identical U applied to H/S/A"
    return WannierOperatorBundle(np.array(b.translations),np.array(b.lattice),
           h,s,p,centers,meta,dict(b.units)).validate()

def synthetic_soc_bundle():
    """Orthonormal spin-1/2 x two-orbital SOC toy; NOT material data."""
    x=np.array([[0,1],[1,0]],complex)
    y=np.array([[0,-1j],[1j,0]],complex)
    z=np.diag([1.,-1.]).astype(complex)
    I=np.eye(2,dtype=complex)
    I4=np.eye(4,dtype=complex)
    kron=np.kron
    h0=.65*kron(z,I)+.41*kron(x,I)+.22*kron(y,x)+.105*kron(I,z)
    plus=.13*kron(z,I)+.075j*kron(y,z)
    h_r=np.stack([plus.conj().T,h0,plus])
    spin_r=np.zeros((3,3,4,4),complex)
    spin_r[1]=np.stack([kron(I,v) for v in (x,y,z)])
    p0=np.stack([np.diag([.12,.12,.63,.63]),.21*I4,.36*I4]).astype(complex)
    p0[0]+=.035*kron(x,I)
    pp=np.zeros((3,4,4),complex)
    pp[0]=.012*kron(z,I)
    pp[1]=.009j*kron(y,x)
    position_r=np.stack([pp.conj().transpose(0,2,1),p0,pp])
    centers=np.diagonal(p0,axis1=-2,axis2=-1).T.real
    meta={"basis_ids":{key:"toy-spin-orbital-WF-v1" for key in ("hamiltonian","spin","position")},
          "spin_kind":"projected_physical_pauli",
          "spin_evidence":"synthetic explicit physical Pauli spin x orbital model",
          "orthonormal_basis":True,"position_convention":"full_connection_including_centers"}
    return WannierOperatorBundle(np.array([[-1,0],[0,0],[1,0]]),
             np.array([[3.4,0,0],[0,4.1,0]]),h_r,spin_r,position_r,centers,meta,dict(UNITS)).validate()
