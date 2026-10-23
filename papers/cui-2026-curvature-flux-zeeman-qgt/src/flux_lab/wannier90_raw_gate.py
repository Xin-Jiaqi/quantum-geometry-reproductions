"""Round33: strict *native* Wannier90 dense hr.dat/r.dat parser with W-S weights.

Official output conventions:
hr.dat: comment, num_wann, nrpts, nrpts Wigner-Seitz degeneracy
integers (15/line normally), then nrpts*num_wann**2 records
Rx Ry Rz m(1-based) n(1-based) ReH ImH.
r.dat: comment, num_wann, nrpts, NO degeneracy header,
and records Rx Ry Rz m n ReX ImX ReY ImY ReZ ImZ.
The matching hr.dat degeneracies supply the r.dat weights.
Under write_ndegen_applied=true both are already weighted and hr degens
are all 1. W-S distance cell-dependent wsvec shifts cannot be silently
omitted; this importer rejects those ambiguous data.
No true source files or Bloch -> Wannier physical-spin transformation
are supplied. These parsers DO NOT certify real material KME.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import hashlib
import numpy as np

def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def _read_lines(path, max_bytes=250_000_000):
    p=Path(path)
    if not p.is_file() or p.stat().st_size>max_bytes:
        raise ValueError("source missing or exceeds bounded maximum file size")
    try:lines=p.read_text(encoding="utf-8").splitlines()
    except UnicodeError as e:raise ValueError("native dense Wannier90 text required") from e
    if len(lines)<4:raise ValueError("Wannier90 file too short")
    return lines

def _header(lines):
    try:n=int(lines[1].strip());nr=int(lines[2].strip())
    except (ValueError,IndexError) as e:raise ValueError("invalid Wannier90 count header") from e
    if not 2<=n<=128 or not 1<=nr<=50_000 or n*n*nr>3_000_000:
        raise ValueError("Wannier90 dimensions outside validated bounded range")
    return n,nr

def _integer(s):
    try:
        i=int(s)
    except ValueError as e:raise ValueError("integer Wannier index or translation required") from e
    return i

def _record(line,width):
    fields=line.split()
    if len(fields)!=width:raise ValueError("malformed Wannier matrix entry columns")
    values=tuple(_integer(x) for x in fields[:5])
    try:
        # Fortran D exponent is also common in text exports.
        floats=tuple(float(x.replace("D","E").replace("d","e")) for x in fields[5:])
    except ValueError as e:raise ValueError("invalid Wannier matrix numeric entry") from e
    if not np.all(np.isfinite(floats)):
        raise ValueError("nonfinite Wannier matrix elements")
    return values,floats

def _matrix_records(lines,offset,n,nr,width,ncomp):
    count=nr*n*n
    if len(lines)!=offset+count:
        raise ValueError("incomplete/extra Wannier matrix records")
    rvec=[];data=np.zeros((nr,ncomp,n,n),complex)
    for ir in range(nr):
        seen=set();rr=None
        for j in range(n*n):
            (rx,ry,rz,m,nn),floats=_record(lines[offset+ir*n*n+j],width)
            r=(rx,ry,rz)
            if rz!=0:raise ValueError("two-dimensional source requires Rz=0")
            if rr is None:rr=r
            if r!=rr:raise ValueError("mixed translations within Wannier R block")
            if not (1<=m<=n and 1<=nn<=n):
                raise ValueError("one-based Wannier orbital index out of range")
            if (m,nn) in seen:raise ValueError("duplicate matrix entry")
            seen.add((m,nn))
            for a in range(ncomp):
                data[ir,a,m-1,nn-1]=complex(floats[2*a],floats[2*a+1])
        if len(seen)!=n*n:raise ValueError("incomplete orbital matrix block")
        if rr in rvec:raise ValueError("duplicate R block")
        rvec.append(rr)
    return np.asarray(rvec,int),data

def parse_hr(path):
    """Exact dense 2D Wannier90_hr.dat syntax, not TB in sparse format."""
    lines=_read_lines(path);n,nr=_header(lines)
    deg=[];pos=3
    while len(deg)<nr:
        if pos>=len(lines):raise ValueError("truncated Wannier degeneracy block")
        tokens=lines[pos].split()
        if not tokens:raise ValueError("empty Wannier degeneracy row")
        deg.extend(_integer(t) for t in tokens);pos+=1
        if len(deg)>nr:raise ValueError("degeneracy length exceeds nrpts")
    d=np.asarray(deg,int)
    if np.any(d<1):raise ValueError("Wigner-Seitz degeneracy must be positive")
    R,mats=_matrix_records(lines,pos,n,nr,7,1)
    return dict(R=R,raw=mats[:,0],degeneracy=d,nwann=n,nrpts=nr,
                source_sha256=_sha(path))

def parse_r(path):
    """Native r.dat without degeneracy block, in Cartesian Angstrom declared."""
    lines=_read_lines(path);n,nr=_header(lines)
    R,mats=_matrix_records(lines,3,n,nr,11,3)
    return dict(R=R,raw=mats,nwann=n,nrpts=nr,
                source_sha256=_sha(path))

@dataclass(frozen=True)
class RawWannierGate:
    R:np.ndarray
    degeneracy:np.ndarray
    h_r:np.ndarray
    position_r:np.ndarray
    centers:np.ndarray
    h_sha256:str
    r_sha256:str
    provenance_status:str="PARSED_H_R_SPIN_NOT_AUTHENTICATED"
    physical_spin_available:bool=False
    original_paper_status:str="PARTIAL"
    kme_ready:bool=False

    def fourier(self,q):
        q=np.asarray(q,float)
        if q.shape!=(2,) or not np.all(np.isfinite(q)):
            raise ValueError("fractional 2D k point required")
        phase=np.exp(2j*np.pi*(self.R[:,:2]@q))
        H=np.einsum("r,rij->ij",phase,self.h_r,optimize=True)
        A=np.einsum("r,raij->aij",phase,self.position_r,optimize=True)
        return H,A

def staged_h_r_import(hr_path,r_path,*,use_ws_distance,write_ndegen_applied,
                      wsvec_present=False,transl_inv_full=False):
    """Reads actual *bytes* but blocks spin/source and incompatible W-S flags.

    Flags must be independently established from original .win/.wout;
    declaring them is not independently authenticating physical settings.
    """
    flags=(use_ws_distance,write_ndegen_applied,wsvec_present,transl_inv_full)
    if not all(isinstance(v,bool) for v in flags):
        raise ValueError("explicit boolean Wannier90 source flags required")
    if use_ws_distance and not write_ndegen_applied:
        raise ValueError("unexpanded ws_distance requires wsvec mapping; not supported")
    if wsvec_present and not write_ndegen_applied:
        raise ValueError("wsvec-required source not supported")
    if transl_inv_full and use_ws_distance and not write_ndegen_applied:
        raise ValueError("full translationally invariant position needs expanded W-S")
    h=parse_hr(hr_path);r=parse_r(r_path)
    if h["nwann"]!=r["nwann"] or h["nrpts"]!=r["nrpts"] or not np.array_equal(h["R"],r["R"]):
        raise ValueError("hr.dat/r.dat translations and orbital order mismatch")
    if write_ndegen_applied and np.any(h["degeneracy"]!=1):
        raise ValueError("write_ndegen_applied requires all hr degeneracies one")
    R=h["R"];idx={tuple(v):i for i,v in enumerate(R)}
    if (0,0,0) not in idx:raise ValueError("onsite R=0 missing")
    for i,v in enumerate(R):
        j=idx.get(tuple(-v))
        if j is None:raise ValueError("missing -R Hermitian partner")
        if h["degeneracy"][i]!=h["degeneracy"][j]:
            raise ValueError("inverse translation degeneracy mismatch")
    w=h["degeneracy"].astype(float)
    # r.dat carries no weights. Both must use the same hr W-S multiplicity.
    H=h["raw"]/w[:,None,None]
    A=r["raw"]/w[:,None,None,None]
    for X in (H,A):
        for i,v in enumerate(R):
            j=idx[tuple(-v)]
            if not np.allclose(X[i],np.swapaxes(X[j].conj(),-1,-2),atol=1e-8,rtol=0):
                raise ValueError("Hermitian -R matrix pairing violated")
    centers=np.stack([np.real(np.diag(A[idx[(0,0,0)],a])) for a in range(3)],axis=1)
    if np.max(abs(np.imag(np.diagonal(A[idx[(0,0,0)]],axis1=-2,axis2=-1))))>1e-9:
        raise ValueError("onsite orbital centers must be real")
    return RawWannierGate(R=R,degeneracy=h["degeneracy"],h_r=H,
       position_r=A,centers=centers,h_sha256=h["source_sha256"],
       r_sha256=r["source_sha256"])
