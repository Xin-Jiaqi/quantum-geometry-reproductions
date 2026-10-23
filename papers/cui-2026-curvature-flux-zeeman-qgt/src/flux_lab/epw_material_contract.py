"""Round32: fail-closed normalized EPW scattering input and material evidence gate.

Not an EPW .epmatwp/.epmat binary reader. Numpy NPZ is an explicitly
separate *normalized export* of matched electronic/phonon scattering input;
manifest evidence must point to actual source files for candidate datasets.
Checksums authenticate bytes, NOT physical provenance, spin, DFPT or Wannier
gauge. No dataset is certified as ab-initio by this code.

Convention: g[k,q,nu,m_final,n_initial] in eV with phase-fixed reverse
phonon modes; each phonon branch paired at +/-q. h[k] Hermitian eV;
S[k,cart] physical projected Pauli; A[k,cart] Angstrom including centers.
A k-dependent unitary changes the connection inhomogeneously, so only
the H/g spectral rate is tested for independent k-dependent rotations.
Finite Gaussian energy broadening does not obey *exact* on-shell detailed
balance. This is NOT a conserving scattering vertex or material KME.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import numpy as np

SCHEMA="cui-epw-normalized-scattering-v1"
UNITS={"electronic_energy":"eV","phonon_energy":"eV","coupling":"eV",
       "position_connection":"angstrom","physical_spin":"dimensionless_pauli",
       "reciprocal":"fractional"}
ROLES={"electronic_scf","phonon_dfpt","electron_phonon_complex",
       "wannier_hamiltonian","wannier_physical_spin","wannier_position",
       "exporter_mapping"}
HBAR_EV_S=6.582119569e-16

def _hash(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for piece in iter(lambda:f.read(1<<18),b""):h.update(piece)
    return h.hexdigest()

def _within(base,relative):
    if not isinstance(relative,str) or not relative or Path(relative).is_absolute():
        raise ValueError("source paths must be relative and inside manifest directory")
    p=(base/relative).resolve()
    if not p.is_relative_to(base.resolve()):
        raise ValueError("source path traversal or symlink outside manifest directory")
    if not p.is_file():
        raise ValueError("evidence source file absent: "+relative)
    return p

def _finite(name,x,shape=None):
    x=np.asarray(x)
    if shape is not None and x.shape!=shape:
        raise ValueError(name+" shape mismatch")
    if not np.all(np.isfinite(x)):raise ValueError(name+" contains nonfinite entries")
    return x

def _coords(x,name):
    x=_finite(name,x)
    if x.ndim!=2 or x.shape[1]!=2 or len(x)==0:
        raise ValueError(name+" must have (N,2) fractional coordinates")
    if np.max(abs(x-np.mod(x,1)))>1e-9:
        raise ValueError(name+" must be fractional in [0,1)")
    if any(np.max(abs(a-b))<1e-8 for i,a in enumerate(x) for b in x[i+1:]):
        raise ValueError("duplicate "+name+" vector")
    return x

def _mesh_maps(k,q):
    """Explicit full closure of k+q and q->-q, no nearest-neighbor guess."""
    plus=np.zeros((len(k),len(q)),int)
    inv=np.zeros(len(q),int)
    for iq,v in enumerate(q):
        qminus=np.mod(-v,1.)
        matches=[j for j,w in enumerate(q) if np.max(abs(np.mod(w-qminus+.5,1.)-.5))<1e-8]
        if len(matches)!=1:raise ValueError("inverse phonon -q is missing/ambiguous")
        inv[iq]=matches[0]
        for ik,u in enumerate(k):
            target=np.mod(u+v,1.)
            candidates=[j for j,w in enumerate(k) if np.max(abs(np.mod(w-target+.5,1.)-.5))<1e-8]
            if len(candidates)!=1:raise ValueError("k+q not represented on electronic mesh")
            plus[ik,iq]=candidates[0]
    return plus,inv

def validate_normalized(arrays,meta,*,folder=None):
    if not isinstance(meta,dict) or meta.get("schema")!=SCHEMA:
        raise ValueError("unsupported explicit scattering schema")
    if meta.get("stage") not in ("synthetic_fixture","ab_initio_candidate"):
        raise ValueError("must distinguish synthetic vs ab_initio_candidate")
    if meta.get("units")!=UNITS:raise ValueError("explicit compatible units required")
    if not isinstance(meta.get("material_id"),str) or not meta["material_id"].strip():
        raise ValueError("material_id required")
    b=meta.get("basis")
    if (not isinstance(b,dict) or set(b.get("gauge_ids",{}))!=
        {"H","g","physical_spin","position"} or
        len(set(b["gauge_ids"].values()))!=1 or
        not all(isinstance(v,str) and v.strip() for v in b["gauge_ids"].values())):
        raise ValueError("shared electronic Wannier gauge IDs required")
    if b.get("spin_kind")!="projected_physical_pauli":
        raise ValueError("not an orbital pseudospin substitute")
    if b.get("position_convention")!="full_connection_including_centers":
        raise ValueError("connection must include centers exactly once")
    if b.get("phonon_phase")!="paired_q_inverse_mode":
        raise ValueError("phonon +/-q branch phase pairing required for conjugation check")
    if b.get("g_indices")!="k,q,nu,m_final,n_initial":
        raise ValueError("electron-phonon g final/input index convention required")
    if meta.get("polar_correction") not in ("not_applicable","included_long_range","short_long_separated"):
        raise ValueError("polar long-range treatment must be declared, not assumed")
    required={"k_frac","q_frac","lattice_angstrom","h_eV","g_eV","phonon_eV",
              "spin_pauli","connection_angstrom","phonon_active"}
    if set(arrays)!=required:raise ValueError("normalized array set incomplete or unexpected")
    k=_coords(arrays["k_frac"],"k_frac")
    q=_coords(arrays["q_frac"],"q_frac")
    lattice=_finite("lattice_angstrom",arrays["lattice_angstrom"],(2,3))
    if np.linalg.matrix_rank(lattice)<2:raise ValueError("degenerate direct lattice")
    h=_finite("h_eV",arrays["h_eV"])
    if h.ndim!=3 or h.shape[0]!=len(k) or h.shape[1]!=h.shape[2] or h.shape[-1]<2:
        raise ValueError("h_eV must be (Nk,Nband,Nband)")
    nb=h.shape[-1];nk=len(k);nq=len(q)
    g=_finite("g_eV",arrays["g_eV"])
    if g.ndim!=5 or g.shape[:2]!=(nk,nq) or g.shape[3:]!=(nb,nb) or g.shape[2]<1:
        raise ValueError("g_eV must be (Nk,Nq,Nmode,Nband,Nband)")
    nm=g.shape[2]
    w=_finite("phonon_eV",arrays["phonon_eV"],(nq,nm)).astype(float)
    active=np.asarray(arrays["phonon_active"])
    if active.shape!=(nq,nm) or active.dtype!=np.dtype(bool):
        raise ValueError("phonon_active must be explicit Boolean mask")
    if np.any(w<0) or np.any(active&(w<=0)):
        raise ValueError("imaginary or nonpositive active phonon energy")
    if np.any((~active)&(w!=0)):
        raise ValueError("inactive mode must be explicit zero-frequency excluded channel")
    s=_finite("spin_pauli",arrays["spin_pauli"],(nk,3,nb,nb))
    a=_finite("connection_angstrom",arrays["connection_angstrom"],(nk,2,nb,nb))
    for name,z in (("h",h),("spin",s),("position",a)):
        if np.max(abs(z-np.swapaxes(z.conj(),-1,-2)))>1e-9:
            raise ValueError(name+" not Hermitian in supplied electronic basis")
    for operator in s.reshape((-1,nb,nb)):
        if np.linalg.norm(operator,2)>1+1e-9:
            raise ValueError("physical projected Pauli operator norm exceeds one")
    plus,inv=_mesh_maps(k,q)
    for iq in range(nq):
        if np.max(abs(w[iq]-w[inv[iq]]))>1e-9 or np.any(active[iq]!=active[inv[iq]]):
            raise ValueError("paired phonon +/-q branch frequencies/masks differ")
        for ik in range(nk):
            j=plus[ik,iq]
            # Phase fixed under +/-q normalized export contract.
            if np.max(abs(g[ik,iq]-np.swapaxes(g[j,inv[iq]].conj(),-1,-2)))>1e-8:
                raise ValueError("electron-phonon g missing inverse q Hermiticity/paired phase")
        if np.any(~active[iq]) and np.max(abs(g[:,iq,~active[iq]]))>1e-12:
            raise ValueError("inactive acoustic mode must have zero exported g")
    evidence=meta.get("evidence",[])
    if not isinstance(evidence,list):raise ValueError("evidence list required")
    if meta["stage"]=="ab_initio_candidate":
        if folder is None:raise ValueError("real candidate requires readable source provenance")
        roles=set()
        for row in evidence:
            if not isinstance(row,dict) or set(row)!={"role","path","sha256"}:
                raise ValueError("source evidence requires role,path,sha256")
            if row["role"] not in ROLES or row["role"] in roles:
                raise ValueError("unknown/duplicate raw-source role")
            roles.add(row["role"])
            if not isinstance(row["sha256"],str) or len(row["sha256"])!=64:
                raise ValueError("source SHA256 required")
            if _hash(_within(Path(folder),row["path"]))!=row["sha256"]:
                raise ValueError("source evidence SHA256 mismatch")
        if roles!=ROLES:raise ValueError("complete original source roles required for candidate")
        if not isinstance(meta.get("exporter_version"),str) or not meta["exporter_version"].strip():
            raise ValueError("mapped EPW exporter provenance/version required")
        status="CANDIDATE_BYTES_VALIDATED_PHYSICAL_AUTHENTICITY_UNVERIFIED"
    else:
        if evidence:raise ValueError("synthetic fixture must not counterfeit external source evidence")
        status="SYNTHETIC_FIXTURE_ONLY"
    return dict(status=status,nk=nk,nq=nq,nband=nb,nmode=nm,
        kplusq=plus,q_inverse=inv,active_mode_count=int(active.sum()),
        real_material_verified=False,paper_status="PARTIAL")

def load_bundle(directory):
    """Explicit NPZ+JSON normalized export, not raw EPW reader."""
    folder=Path(directory).resolve()
    meta_file=folder/"manifest.json";array_file=folder/"normalized.npz"
    if not meta_file.is_file() or not array_file.is_file():
        raise ValueError("missing normalized.npz or manifest.json")
    meta=json.loads(meta_file.read_text(encoding="utf-8"))
    digest=meta.get("normalized_sha256")
    if not isinstance(digest,str) or _hash(array_file)!=digest:
        raise ValueError("normalized NPZ SHA256 mismatch")
    with np.load(array_file,allow_pickle=False) as bundle:
        data={k:bundle[k] for k in bundle.files}
    verdict=validate_normalized(data,meta,folder=folder)
    return data,meta,verdict

def scatter_rates(arrays,verdict,*,beta=40.0,mu=0.0,sigma_eV=.08):
    """Pauli-blocked finite-width golden-rule *diagnostic*, not SI transport.

    Rates 1/s of each band scattering event; Gaussian delta in 1/eV.
    Finite sigma is an off-shell approximation and must NOT be declared
    exact detailed balance or microscopic self-energy/vertex.
    """
    if not np.isfinite(beta) or beta<=0 or not np.isfinite(mu):
        raise ValueError("positive beta (1/eV) and finite mu required")
    if not np.isfinite(sigma_eV) or sigma_eV<=0:
        raise ValueError("positive energy broadening in eV required")
    kplusq=verdict["kplusq"];h=arrays["h_eV"];g=arrays["g_eV"]
    w=arrays["phonon_eV"];mask=arrays["phonon_active"]
    nk,nq,nmode,nb,_=g.shape
    eigs=[];eigvecs=[]
    for a in h:
        e,U=np.linalg.eigh(a);eigs.append(e);eigvecs.append(U)
    e=np.asarray(eigs);u=np.asarray(eigvecs)
    f=1/(1+np.exp(np.clip(beta*(e-mu),-700,700)))
    bare=np.zeros((2,nk,nq,nmode,nb,nb),float)
    for ik in range(nk):
        for iq in range(nq):
            jk=kplusq[ik,iq]
            for nu in range(nmode):
                if not mask[iq,nu]:continue
                ww=w[iq,nu]
                gband=u[jk].conj().T@g[ik,iq,nu]@u[ik]
                delta=e[jk,:,None]-e[ik,None,:]
                nB=1/np.expm1(beta*ww)
                gaussian=lambda x:np.exp(-.5*(x/sigma_eV)**2)/(sigma_eV*np.sqrt(2*np.pi))
                pref=(2*np.pi/HBAR_EV_S)*abs(gband)**2
                bare[0,ik,iq,nu]=pref*nB*gaussian(delta-ww) # absorption
                bare[1,ik,iq,nu]=pref*(nB+1)*gaussian(delta+ww) # emission
    blocked=np.zeros_like(bare)
    for ik in range(nk):
        for iq in range(nq):
            blocked[:,ik,iq]=bare[:,ik,iq]*(1-f[kplusq[ik,iq]])[None,None,:,None]
    return dict(bare=bare,blocked=blocked,eigenvalues_eV=e,fermi=f,
                unit="s^-1",finite_broadening_exact_KMS=False,
                real_material_verified=False)

def make_synthetic_fixture():
    """A three-k, two spin-band fixture with exact reciprocal-mode pairing."""
    nk=3;nq=3;nb=2;nmode=1
    k=np.array([[j/3,0.] for j in range(3)])
    q=k.copy()
    pauli_x=np.array([[0,1],[1,0]],complex)
    pauli_y=np.array([[0,-1j],[1j,0]],complex)
    pauli_z=np.diag([1.,-1.]).astype(complex)
    h=np.array([np.diag([-.15+.03*j,.18+.02*j]) for j in range(3)],complex)
    g=np.zeros((nk,nq,nmode,nb,nb),complex)
    # q=0 Hermitian; q=1/3 and 2/3 paired with exact reversed endpoints.
    for j in range(3):g[j,0,0]=.015*pauli_x+.011*np.eye(2)
    for j in range(3):
        M=(.02+.002j*j)*pauli_x+.006j*pauli_y+.005*np.eye(2)
        g[j,1,0]=M
        g[(j+1)%3,2,0]=M.conj().T
    s=np.tile(np.stack([pauli_x,pauli_y,pauli_z])[None],(3,1,1,1))
    a=np.zeros((3,2,nb,nb),complex)
    a[:,0]=np.diag([.08,.22])
    a[:,1]=np.diag([.02,.12])
    arrays=dict(k_frac=k,q_frac=q,lattice_angstrom=np.array([[3.,0.,0.],[0.,3.,0.]]),
        h_eV=h,g_eV=g,phonon_eV=np.array([[.07],[.07],[.07]]),
        spin_pauli=s,connection_angstrom=a,phonon_active=np.ones((3,1),bool))
    meta=dict(schema=SCHEMA,stage="synthetic_fixture",material_id="SYNTHETIC_NOT_A_MATERIAL",
        units=UNITS,basis=dict(gauge_ids=dict.fromkeys(("H","g","physical_spin","position"),"SYN_R32"),
            spin_kind="projected_physical_pauli",
            position_convention="full_connection_including_centers",
            phonon_phase="paired_q_inverse_mode",
            g_indices="k,q,nu,m_final,n_initial"),
        polar_correction="not_applicable",evidence=[])
    return arrays,meta

def rotate_electronic_gauge(arrays,unitaries,plus):
    """Test k-dependent band H/g covariance, *not* A's derivative connection."""
    u=np.asarray(unitaries,complex);nb=arrays["h_eV"].shape[-1]
    if u.shape!=(len(arrays["k_frac"]),nb,nb) or any(
        not np.allclose(M.conj().T@M,np.eye(nb),atol=1e-11) for M in u):
        raise ValueError("unitary matrices required at each k")
    out={name:np.array(value,copy=True) for name,value in arrays.items()}
    for ik in range(len(u)):
        out["h_eV"][ik]=u[ik].conj().T@arrays["h_eV"][ik]@u[ik]
        for iq in range(len(arrays["q_frac"])):
            for nu in range(arrays["g_eV"].shape[2]):
                out["g_eV"][ik,iq,nu]=(u[plus[ik,iq]].conj().T@
                                         arrays["g_eV"][ik,iq,nu]@u[ik])
    # A is NOT updated; claiming an independent-k connection gauge would be wrong.
    return out
