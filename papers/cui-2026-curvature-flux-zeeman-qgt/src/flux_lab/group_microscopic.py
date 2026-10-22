"""Round16 microscopic bridge: full 31 2D MPG local interband symmetry.

The pinned magnetic_symmetry catalogue records REAL-SPACE operations R.
Momentum is transformed by Q=(-1)^anti R. The physical axial spin
transforms by L=(-1)^anti det(R) diag(R,1). For A=U K^anti:
chi_JB(Qk,z)=Q chi_JB(k,(-1)^anti z) L.T;
alpha_ME(Qk,z)=(-1)^anti L alpha_ME(k,(-1)^anti z) Q.T.
At eta>0 the -z branch is advanced, not retarded(-omega+i eta).

All group checks are for synthetic four-band nondegenerate models.
C3/C6 local algebra does NOT supply a properly periodic square-grid BZ.
No collision, real-material magnetic response, Kramers doublet, or SI claim.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .core import PAULI
from .magnetic_symmetry import load_groups,matrices,SECTORS,allowed
from .causal_response import interband_kubo,interband_liouville

I=np.eye(2,dtype=complex)
S=np.stack([np.kron(I,s) for s in PAULI])

@dataclass(frozen=True)
class Operation:
    R:np.ndarray
    anti:bool
    U:np.ndarray

def _generator(name,anti):
    if name=="E":u=I
    elif name=="Mx":u=1j*PAULI[0]
    elif name.startswith("C"):
        theta=2*np.pi/int(name[1:])
        u=np.cos(theta/2)*I-1j*np.sin(theta/2)*PAULI[2]
    else:raise ValueError("unknown generator")
    if anti:u=u@(1j*PAULI[1])
    return Operation(matrices()[name],bool(anti),u)

def momentum_map(g):
    return (-1 if g.anti else 1)*g.R

def spin_axial(g):
    R3=np.eye(3);R3[:2,:2]=g.R
    return ((-1 if g.anti else 1)*np.linalg.det(g.R))*R3

def build_group(name):
    groups=load_groups()
    if name not in groups:raise ValueError("not a pinned 2D MPG")
    gens=[_generator(a,t) for a,t in groups[name]["generators"]]
    ops=[Operation(np.eye(2),False,I)]
    def same(a,b):return a.anti==b.anti and np.max(abs(a.R-b.R))<1e-10
    while True:
        new=list(ops)
        for a in ops:
            for g in gens:
                b=Operation(a.R@g.R,bool(a.anti^g.anti),
                            a.U@(g.U.conj() if a.anti else g.U))
                found=next((p for p in new if same(p,b)),None)
                if found is None:new.append(b)
                else:
                    phase=np.trace(found.U.conj().T@b.U)/2
                    if not np.allclose(b.U,phase*found.U,atol=2e-10):
                        raise AssertionError("spinor representation not projectively closed")
        if len(new)==len(ops):break
        if len(new)>40:raise ValueError("nonfinite point group")
        ops=new
    if len(ops)!=len(groups[name]["operations"]):
        raise AssertionError("wrong point-group order")
    for g in ops:
        got=np.stack([g.U@(s.conj() if g.anti else s)@g.U.conj().T for s in PAULI])
        want=np.einsum("bc,cij->bij",spin_axial(g).T,PAULI)
        if not np.allclose(got,want,atol=2e-10):
            raise AssertionError("physical spin and MPG parity disagree")
    return ops

def _seed(k):
    k=np.asarray(k,float)
    x,y,z=PAULI
    terms=[
        (2.6*np.kron(z,I),"constant",(0,0)),
        (.43*np.kron(I,x),"sin",(1,0)),
        (.27*np.kron(I,x),"sin",(0,1)),
        (-.27*np.kron(I,y),"sin",(1,0)),
        (.31*np.kron(x,y),"cos",(0,1)),
        (.35*np.kron(z,z),"sin",(1,1)),
        (.21*np.kron(y,x),"cos",(1,-1)),
        (.44*np.kron(I,y),"sin",(0,1)),
        (.29*np.kron(x,z),"cos",(1,0)),
        (.14*np.kron(I,I),"cos",(1,1)),
        (.39*np.kron(y,z),"sin",(1,-1))]
    H=np.zeros((4,4),complex);V=np.zeros((2,4,4),complex)
    for A,typ,q in terms:
        if typ=="constant":H+=A;continue
        q=np.asarray(q,float);phase=q@k
        f,fp=(np.sin(phase),np.cos(phase)) if typ=="sin" else (np.cos(phase),-np.sin(phase))
        H+=f*A;V+=q[:,None,None]*fp*A
    return H,V

def symmetrized_model(k,group):
    """Exact group twirl of a generic periodic seed and its analytic gradient."""
    ops=build_group(group) if isinstance(group,str) else group
    H=np.zeros((4,4),complex);V=np.zeros((2,4,4),complex)
    for g in ops:
        Q=momentum_map(g)
        h,v=_seed(Q.T@np.asarray(k))
        U4=np.kron(I,g.U)
        trans=lambda A:U4@(A.conj() if g.anti else A)@U4.conj().T
        H+=trans(h)
        V+=np.einsum("ab,bij->aij",Q,np.stack([trans(w) for w in v]))
    return H/len(ops),V/len(ops),S.copy()

def local_group_residual(name,k,z=.051+.011j,mu=.17,temperature=.17,liouville=False):
    """Independent pair response transformed for EVERY operation of this MPG."""
    ops=build_group(name)
    h,v,s=symmetrized_model(k,ops)
    method=interband_liouville if liouville else interband_kubo
    errors=np.zeros(4)
    for g in ops:
        Q=momentum_map(g);L=spin_axial(g)
        hp,vp,sp=symmetrized_model(Q@np.asarray(k),ops)
        U4=np.kron(I,g.U)
        trans=lambda A:U4@(A.conj() if g.anti else A)@U4.conj().T
        errors[0]=max(errors[0],np.max(abs(hp-trans(h))))
        predv=np.einsum("ab,bij->aij",Q,np.stack([trans(w) for w in v]))
        errors[1]=max(errors[1],np.max(abs(vp-predv)))
        chi,alpha=method(h,v,s,z=-z if g.anti else z,mu=mu,temperature=temperature)
        chi1,alpha1=method(hp,vp,sp,z=z,mu=mu,temperature=temperature)
        errors[2]=max(errors[2],np.max(abs(chi1-Q@chi@L.T)))
        errors[3]=max(errors[3],np.max(abs(alpha1-(-1 if g.anti else 1)*L@alpha@Q.T)))
    return dict(H=float(errors[0]),velocity=float(errors[1]),
                chi=float(errors[2]),alpha=float(errors[3]),order=len(ops))

def sector_nullspace(name,parity,component="inplane"):
    """Exact 2x3 JB even/odd group constraint, no response fitting."""
    if parity not in ("even","odd"):raise ValueError("even or odd required")
    if component not in ("inplane","full"):raise ValueError("inplane or full required")
    bs=[]
    for a in range(2):
        for b in range(2 if component=="inplane" else 3):
            x=np.zeros((2,3));x[a,b]=1;bs.append(x)
    constraints=[]
    for g in build_group(name):
        factor=-1 if parity=="odd" and g.anti else 1
        constraints.append(np.stack([(x-factor*momentum_map(g)@x@spin_axial(g).T).ravel()
                                     for x in bs],axis=1))
    values=np.linalg.svd(np.vstack(constraints),compute_uv=False)
    return len(bs)-int(np.sum(values>1e-9))

def verify_catalog_projectors():
    out=[]
    for name,group in load_groups().items():
        for parity,eta in (("even",1),("odd",-1)):
            local=sector_nullspace(name,parity)
            jahns=sum(allowed(group,s) for s,(_,t) in SECTORS.items() if t==eta)
            out.append(dict(group=name,parity=parity,microscopic_rank=local,
                            jahnscreen_rank=jahns,agree=local==jahns))
    return out

def project_frequency_tensor(name,X,parity):
    if parity not in ("even","odd"):raise ValueError("even or odd required")
    X=np.asarray(X,complex)
    if X.shape!=(2,3):raise ValueError("2x3 J/B tensor required")
    ops=build_group(name)
    return sum(((-1 if parity=="odd" and g.anti else 1)*
                momentum_map(g)@X@spin_axial(g).T for g in ops),
               np.zeros((2,3),complex))/len(ops)

def clean_bz_even_odd(name,n=8,z=.043,mu=.17,temperature=.17):
    """Square-torus k integration; C3/C6 require a hexagonal BZ and are rejected."""
    if n<4 or n%2:raise ValueError("even mesh >=4 required")
    group=load_groups()[name]
    if any(token not in ("E","Mx","C2","C4") for token,_ in group["generators"]):
        raise ValueError("hexagonal reciprocal BZ required for C3/C6")
    ops=build_group(name)
    ks=-np.pi+(np.arange(n)+.5)*2*np.pi/n
    pos=np.zeros((2,3),complex);neg=pos.copy()
    for x in ks:
        for y in ks:
            h,v,s=symmetrized_model([x,y],ops)
            pos+=interband_kubo(h,v,s,z=z,mu=mu,temperature=temperature)[0]
            neg+=interband_kubo(h,v,s,z=-z,mu=mu,temperature=temperature)[0]
    return (pos+neg)/(2*n*n),(pos-neg)/(2*n*n)
