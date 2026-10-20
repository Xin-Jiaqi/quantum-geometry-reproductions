"""IGMC MPG tensor audit using pinned generators from author's group-theory toolkit.

This is an independent Neumann null-space solver using the toolkit's *operations*,
not a copy of its invariant solver. No imported toolkit is needed offline.
Optional live cross-check calls the public toolkit magnetic_tensor_basis API.
The input snapshot is intentionally limited to 122 MPG generators, not the full repo.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
GENERATOR_FILE = ROOT/'data'/'toolkit_mpg_generators_v018.txt'
EXPECTED_ORDERS = tuple(map(int, ('1,2,2,4,2,2,4,2,2,4,2,4,8,4,4,4,4,8,4,4,8,4,4,8,16,8,8,8,4,8,4,4,8,4,8,16,8,8,8,8,16,8,8,8,16,8,8,8,16,8,8,8,16,32,16,16,16,16,16,3,6,6,12,6,6,12,6,6,12,6,12,24,12,12,12,6,12,6,6,12,6,12,24,12,12,12,12,24,12,12,12,24,12,12,12,24,12,12,12,24,48,24,24,24,24,24,12,24,24,48,24,24,48,24,24,48,24,48,96,48,48,48').split(',')))

@dataclass(frozen=True)
class MagneticGroup:
    number:int
    symbol:str
    family:str
    generators:tuple
    order:int


def read_toolkit_snapshot(path:Path=GENERATOR_FILE):
    maps={};groups=[]
    for raw in Path(path).read_text(encoding='utf-8').splitlines():
        if raw.startswith('MAP '):
            fam, rest=raw[4:].split(' ',1)
            mapping={}
            for token in rest.split(';'):
                code,op,flat=token.split(':',2)
                matrix=np.array([float(y) for y in flat.split(',')]).reshape(3,3)
                if not np.allclose(matrix@matrix.T,np.eye(3),atol=2e-7):
                    raise ValueError(f'nonorthogonal toolkit matrix {fam}/{op}')
                mapping[code]=(op,matrix)
            maps[fam]=mapping
        if raw.startswith('RECORDS '):
            for entry in raw[8:].split('|'):
                head,keys=entry.split(':',1)
                idx,sym=head.split('=',1); idx=int(idx)
                family=('cubic_Oh' if idx<=5 or idx>=107 else
                        'tetragonal_D4h' if idx<=59 else 'hexagonal_D6h')
                gen=[]
                for key in keys.split(',') if keys else ():
                    op,matrix=maps[family][key.replace('!','')]
                    gen.append((op,matrix,key.endswith('!')))
                groups.append(MagneticGroup(idx,sym,family,tuple(gen),EXPECTED_ORDERS[idx-1]))
    if len(groups)!=122 or [g.number for g in groups]!=list(range(1,123)):
        raise ValueError('magnetic-point-group snapshot incomplete or out of order')
    return tuple(groups)


def tensor_matrix(rotation,antiunitary,channel='C',protocol='channel_resolved_symmetry'):
    """Mapping of real polar-current <- axial-magnetic-field coefficients."""
    if channel not in {'C','D'}:raise ValueError(channel)
    if protocol=='channel_resolved_symmetry':
        sign=(-1 if channel=='D' else +1) if antiunitary else 1
    elif protocol=='paper_printed_eq14':
        sign=-1 if antiunitary else +1
    else:raise ValueError(protocol)
    r=np.asarray(rotation)
    det=int(round(np.linalg.det(r)))
    transform=np.column_stack([
        (sign*det*r@unit.reshape(3,3)@r.T).ravel()
        for unit in np.eye(9)
    ])
    return transform


def tensor_basis(group, channel='C', protocol='channel_resolved_symmetry'):
    """SVD nullspace, 9×dim columns, for all toolkit MPG generators."""
    constraints=[np.eye(9)-tensor_matrix(r,tr,channel,protocol)
                 for name,r,tr in group.generators]
    if not constraints:return np.eye(9)
    _,svh,vh=np.linalg.svd(np.vstack(constraints),full_matrices=True)
    rank=int(np.sum(svh>1e-7))
    return vh[rank:,:].T


def classify_basis(basis):
    if basis.ndim!=2 or basis.shape[0]!=9:raise ValueError('basis shape')
    if basis.shape[1]==0:return 'forbidden'
    diag_idx=(0,4,8)
    diag=bool(np.max(abs(basis[list(diag_idx),:]))>1e-6)
    off=bool(np.max(abs(basis[[1,2,3,5,6,7],:]))>1e-6)
    return 'both' if diag and off else 'longitudinal_only' if diag else 'transverse_only'


def group_invariance_max(group,basis,channel='C',protocol='channel_resolved_symmetry'):
    return max((np.linalg.norm((np.eye(9)-tensor_matrix(r,tr,channel,protocol))@basis)
                for name,r,tr in group.generators),default=0.)


def all_group_operations(group):
    """Enumerate full magnetic group from generators as (orthogonal R, bool T)."""
    gens=[(np.asarray(r),bool(tr)) for _,r,tr in group.generators]
    ident=np.eye(3)
    def key(m,t):return (tuple(np.round(m,5).ravel()),bool(t))
    seen={key(ident,False):(ident,False)}
    pending=[(ident,False)]
    while pending:
        r,t=pending.pop()
        for g,s in gens:
            nr,nt=r@g,t!=s
            k=key(nr,nt)
            if k not in seen:
                seen[k]=(nr,nt);pending.append((nr,nt))
                if len(seen)>96:raise AssertionError(f'group closure exploded: {group.symbol}')
    return tuple(seen.values())


def closure_order(group):
    """Magnetic group closure count with antiunitary Z2 label."""
    return len(all_group_operations(group))


def character_dimension(group,channel='C',protocol='channel_resolved_symmetry'):
    """Independent irrep-character dimension; no SVD or tensor_matrix.

    For polar current J and axial input B, the 9D map representation has
    character det(R)*[tr(R)]**2. Antiunitary real-amplitude channel character
    eta_C=+1 or eta_D=-1 multiplies odd-time-reversal elements. Reynolds
    character average over the COMPLETE magnetic group gives dim Inv(M).
    """
    if channel not in ('C','D'):raise ValueError(channel)
    if protocol=='channel_resolved_symmetry':eta=1 if channel=='C' else -1
    elif protocol=='paper_printed_eq14':eta=-1
    else:raise ValueError(protocol)
    operations=all_group_operations(group)
    if len(operations)!=group.order:
        raise AssertionError(f'bad closure for {group.symbol}')
    chars=[((eta if tr else 1)*round(np.linalg.det(r))*np.trace(r)**2)
           for r,tr in operations]
    average=float(np.mean(chars))
    dim=int(round(average))
    if not (0<=dim<=9) or abs(dim-average)>2e-5:
        raise AssertionError(f'character sum not integer: {group.symbol}: {average}')
    return dim


def toolkit_live_dimensions(groups):
    """Optional live toolkit API integration; imported only on explicit invocation."""
    from group_theory_operations import magnetic_tensor_basis
    failures=[]
    for group in groups:
        for channel in ('C','D'):
            # polar current J is T-odd, Zeeman B is T-odd,
            # time derivative of B is T-even.
            result=magnetic_tensor_basis(group.number,'polar_vector','axial_vector',
                  output_time_parity='odd',input_time_parity=('odd' if channel=='C' else 'even'))
            independent=tensor_basis(group,channel).shape[1]
            if independent!=result.dimension:
                failures.append((group.number,group.symbol,channel,independent,result.dimension))
    return failures


def canonical_component_basis(basis):
    """Make each independent coefficient a chosen Cartesian matrix component.

    Unlike an arbitrary SVD orthonormal basis, this selects lexicographic pivot
    components and sets them to identity, displaying SM-like component equalities.
    The tensor subspace is unchanged and is still checked by Neumann's principle.
    """
    b=np.asarray(basis,dtype=float)
    d=b.shape[1]
    if d==0:return np.zeros((9,0))
    pivots=[]
    for i in range(9):
        if np.linalg.matrix_rank(b[pivots+[i],:],tol=1e-7)>len(pivots):
            pivots.append(i)
            if len(pivots)==d:break
    if len(pivots)!=d:raise AssertionError('cannot select independent tensor components')
    out=b@np.linalg.inv(b[pivots,:])
    out[np.abs(out)<1e-8]=0.
    if not np.allclose(out[pivots,:],np.eye(d),atol=1e-6):
        raise AssertionError('canonical tensor components not normalized')
    return out
