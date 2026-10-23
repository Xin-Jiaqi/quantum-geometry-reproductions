"""Round34: genuine Wannier90 ASCII .spn.fmt, _u.mat, _u_dis.mat and eig parsers.

Physics: native .spn.fmt stores <psi_nk|sigma_a|psi_mk> with m outer,
n<=m inner and Pauli component a fastest. U(k) and Udis(k) are written
column-major. Wannier spin = V^dag S_Bloch V, V=Udis U, or V=U if isolated.

The optional _v.mat independently verifies the composition.
The independent eigenvalue/hr Fourier oracle verifies matching electronic
Hamiltonians on the exact native ordered k mesh; it cannot authenticate
the DFT spin provenance, phase alignment of externally patched files,
or a physical KME material. No raw Fortran-unformatted .spn decoding.
"""
from __future__ import annotations
from pathlib import Path
import hashlib
import numpy as np
from .wannier90_raw_gate import staged_h_r_import

def _read(path,max_bytes=120_000_000):
    p=Path(path)
    if not p.is_file() or p.stat().st_size>max_bytes:raise ValueError("missing/oversize native ASCII input")
    try:s=p.read_text(encoding="utf-8").splitlines()
    except UnicodeError as e:raise ValueError("formatted ASCII required; convert binary .spn with official w90spn2spn.x") from e
    if len(s)<2:raise ValueError("truncated Wannier90 native file")
    return s

def _digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def _float(v):
    try:x=float(v.replace("D","E").replace("d","e"))
    except ValueError as e:raise ValueError("invalid native numeric entry") from e
    if not np.isfinite(x):raise ValueError("nonfinite native numeric entry")
    return x

def _ints(tokens,count):
    try:values=tuple(int(x) for x in tokens)
    except ValueError as e:raise ValueError("invalid integer native header") from e
    if len(values)!=count:raise ValueError("native header count mismatch")
    if any(z<1 for z in values):raise ValueError("native dimensions must be positive")
    return values

def parse_spn_formatted(path):
    """Official header + 3 upper-triangular complex Pauli components per k."""
    lines=_read(path)
    nb,nk=_ints(lines[1].split(),2)
    if nb>128 or nk>25000 or nk*nb*(nb+1)*3//2>7_000_000:
        raise ValueError("native spn dimensions over bounded parser limit")
    rows=[s for s in lines[2:] if s.strip()]
    if len(rows)!=nk*nb*(nb+1)*3//2:
        raise ValueError("formatted spn missing/extra upper-triangle records")
    s=np.zeros((nk,3,nb,nb),complex)
    cursor=0
    max_imag_diag=0.
    for ik in range(nk):
        for m in range(nb):
            for n in range(m+1):
                for a in range(3):
                    words=rows[cursor].split();cursor+=1
                    if len(words)!=2:raise ValueError("spn must have Re Im per Pauli record")
                    value=complex(_float(words[0]),_float(words[1]))
                    if m==n:
                        max_imag_diag=max(max_imag_diag,abs(value.imag))
                        if abs(value.imag)>2e-5:
                            raise ValueError("physical Pauli diagonal has significant imaginary part")
                        value=complex(value.real)
                    s[ik,a,n,m]=value
                    if m!=n:s[ik,a,m,n]=value.conjugate()
    # Spinor Bloch Pauli operators have eigenvalues within +/-1.
    for block in s.reshape((-1,nb,nb)):
        if np.linalg.norm(block,2)>1+1e-5:
            raise ValueError("spin Pauli operator norm exceeds one; inconsistent data")
    return dict(spin_bloch=s,nband=nb,nk=nk,diag_imag_raw_max=max_imag_diag,
                source_sha256=_digest(path),kind="BLOCH_PAULI_NOT_WANNIER")

def _matrix_file(path,role):
    lines=_read(path)
    if role not in ("u","u_dis","v"):raise ValueError("unknown matrix role")
    hdr=_ints(lines[1].split(),3)
    nk,nwann,third=hdr
    nb=nwann if role=="u" else third
    if role=="u" and third!=nwann:
        raise ValueError("_u.mat must be nwann x nwann")
    if nk>25000 or nb>128 or nwann>128 or nk*nb*nwann>3_000_000:
        raise ValueError("matrix dimensions exceed guarded bound")
    entries=[s for s in lines[2:] if s.strip()]
    if len(entries)!=nk*(1+nb*nwann):
        raise ValueError("native u/u_dis/v k blocks are incomplete")
    coords=np.zeros((nk,3))
    matrices=np.zeros((nk,nb,nwann),complex)
    for ik in range(nk):
        at=ik*(1+nb*nwann)
        words=entries[at].split()
        if len(words)!=3:raise ValueError("native matrix k coordinate needs three components")
        coords[ik]=[_float(w) for w in words]
        values=[]
        for j in range(nb*nwann):
            part=entries[at+1+j].split()
            if len(part)!=2:raise ValueError("native matrix elements need Re Im")
            values.append(complex(_float(part[0]),_float(part[1])))
        matrices[ik]=np.asarray(values).reshape((nb,nwann),order="F")
    if len({tuple(np.round(x%1,8)) for x in coords})!=nk:
        raise ValueError("duplicate native k coordinates")
    eye=np.eye(nwann)
    for mat in matrices:
        if not np.allclose(mat.conj().T@mat,eye,atol=2e-7,rtol=0):
            raise ValueError("native U/V columns not orthonormal")
    return dict(matrix=matrices,k_frac=coords,nk=nk,nwann=nwann,
                nband=nb,source_sha256=_digest(path),role=role)

def parse_u(path):
    return _matrix_file(path,"u")

def parse_u_dis(path):
    return _matrix_file(path,"u_dis")

def parse_v(path):
    return _matrix_file(path,"v")

def parse_eig(path,*,nk,nband):
    lines=_read(path)
    entries=[x for x in lines if x.strip()]
    if len(entries)!=nk*nband:
        raise ValueError("native eig requires exactly nk*nband lines")
    eigs=np.zeros((nk,nband))
    for ik in range(nk):
        for ib in range(nband):
            words=entries[ik*nband+ib].split()
            if len(words)!=3:raise ValueError("native eig expects band k energy")
            first,second=_ints(words[:2],2)
            if first!=ib+1 or second!=ik+1:
                raise ValueError("native eig band/k ordering mismatch")
            eigs[ik,ib]=_float(words[2])
    return dict(energies=eigs,source_sha256=_digest(path))

def reconstruct_spin(spn,u,*,u_dis=None,v=None,eig=None):
    nk=spn["nk"];nb=spn["nband"];nw=u["nwann"]
    if u["nk"]!=nk:raise ValueError("spn and U k counts differ")
    if u_dis is None:
        if nb!=nw:raise ValueError("disentanglement matrix needed for nbands>nwann")
        transform=u["matrix"]
    else:
        if u_dis["nk"]!=nk or u_dis["nband"]!=nb or u_dis["nwann"]!=nw:
            raise ValueError("disentanglement dimensions/spn mismatch")
        if not np.allclose(u["k_frac"],u_dis["k_frac"],atol=2e-8,rtol=0):
            raise ValueError("U and Udis native k order mismatched")
        transform=np.matmul(u_dis["matrix"],u["matrix"])
    if v is not None:
        if (v["nk"]!=nk or v["nband"]!=nb or v["nwann"]!=nw or
            not np.allclose(v["k_frac"],u["k_frac"],atol=2e-8,rtol=0)):
            raise ValueError("independent V dimensions/k order mismatch")
        if not np.allclose(transform,v["matrix"],atol=2e-7,rtol=0):
            raise ValueError("independent V != Udis U; wrong combined gauge")
    spin=np.empty((nk,3,nw,nw),complex)
    for ik in range(nk):
        for a in range(3):
            spin[ik,a]=transform[ik].conj().T@spn["spin_bloch"][ik,a]@transform[ik]
            if not np.allclose(spin[ik,a],spin[ik,a].conj().T,atol=2e-8,rtol=0):
                raise ValueError("projected physical Pauli spin not Hermitian")
            if np.linalg.norm(spin[ik,a],2)>1+2e-5:
                raise ValueError("projected physical spin norm exceeds one")
    h_ref=None
    if eig is not None:
        if eig["energies"].shape!=(nk,nb):raise ValueError("eig and spn Bloch bands/k mismatch")
        h_ref=np.array([transform[i].conj().T@
                        np.diag(eig["energies"][i])@transform[i] for i in range(nk)])
    return dict(spin_wannier_k=spin,H_wannier_k_from_eig=h_ref,
                transform=transform,k_frac=u["k_frac"],nband=nb,nwann=nw,nk=nk,
                physical_pauli_not_orbital_pseudospin=True,
                mathematical_gauge_transform_verified=True,
                physical_source_provenance_authenticated=False,
                material_KME_ready=False,
                original_cui_paper="PARTIAL")

def check_against_native_hr(reconstruction,hr,r,*,flags,tol_eV=2e-4):
    """H_W from Bloch energies versus independent native HR Fourier.

    No spin-real-space or material authenticity certification; only a
    bounded consistency gate on *exact same ordered* Wannier k grid.
    """
    if reconstruction["H_wannier_k_from_eig"] is None:
        raise ValueError("native Bloch eigenvalues required for HR cross-check")
    if not isinstance(flags,dict) or set(flags)!={"use_ws_distance",
            "write_ndegen_applied","wsvec_present","transl_inv_full"}:
        raise ValueError("explicit Round33 flag contract required")
    native=staged_h_r_import(hr,r,**flags)
    if native.h_r.shape[-1]!=reconstruction["nwann"]:
        raise ValueError("Wannier orbital dimension differs between u and hr")
    coords=reconstruction["k_frac"]
    if np.max(abs(coords[:,2]))>1e-8:
        raise ValueError("two-dimensional HR reference requires kz=0")
    comparisons=[]
    for i,k in enumerate(coords):
        h,_=native.fourier(k[:2])
        comparisons.append(float(np.max(abs(h-reconstruction["H_wannier_k_from_eig"][i]))))
    residual=max(comparisons)
    if residual>tol_eV:
        raise ValueError("independently reconstructed HR and Bloch eig/U mismatch")
    return dict(h_from_eig_vs_hr_max_eV=residual,
                parsed_native_H_r=True,
                native_spin_is_authenticated=False,
                physical_spin_source_provenance_authenticated=False,
                material_KME_ready=False,original_paper="PARTIAL")
