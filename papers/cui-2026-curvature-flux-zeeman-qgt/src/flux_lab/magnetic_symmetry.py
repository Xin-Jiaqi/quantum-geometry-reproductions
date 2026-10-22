"""2D magnetic point groups: source Jahn-character screening vs full tensor projection.

Toolkit source: BSD-3-Clause, 551a0daf8d4ba5a5f24b2811694bfe438616de48.
Group data are a derived generator-only snapshot; no external toolkit dependency.

Assumed *static real sector coefficient* action for a current (polar) driven by
an in-plane magnetic field (axial): eta**theta det(R3) R2 X R2.T.
eta=-1 for source Jahn ae (gN/OmegaA), +1 for e (gA/OmegaN).
This is a candidate *Jahn tensor* rule, not proof of dynamic antiunitary
Kubo/Onsager law. The latter must be independently demonstrated.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np

REGISTRY = Path(__file__).resolve().parents[2] / "data" / "toolkit_2d_generators.json"
SECTORS = {"gN": ("symmetric", -1), "gA": ("symmetric", +1),
           "omegaN": ("antisymmetric", +1), "omegaA": ("antisymmetric", -1)}


def matrices():
    return {"E": np.eye(2), "Mx": np.diag([-1., 1.]),
            **{f"C{n}": np.array([[np.cos(2*np.pi/n), -np.sin(2*np.pi/n)],
                                   [np.sin(2*np.pi/n), np.cos(2*np.pi/n)]])
               for n in (2, 3, 4, 6)}}


def closure(generators, tol=1e-9):
    """Close a 2D magnetic group using true Cartesian multiplication and T XOR."""
    m = matrices()
    operations = [(np.eye(2), False)]
    for name, theta in generators:
        if name not in m or theta not in (0, 1):
            raise ValueError('unrecognized generator')
    def same(x, y):
        return x[1] == y[1] and np.linalg.norm(x[0]-y[0]) < tol
    while True:
        nxt = list(operations)
        for a, ta in operations:
            for name, tb in generators:
                x = (a @ m[name], bool(ta ^ bool(tb)))
                if not any(same(x, o) for o in nxt):
                    nxt.append(x)
        if len(nxt)==len(operations):
            return operations
        if len(nxt)>100: raise ValueError('nonfinite group, check embedding')
        operations=nxt


def load_groups():
    data=json.loads(REGISTRY.read_text(encoding='utf-8'))
    if data['schema'] != '2d-mpg-generators-v1': raise ValueError('registry schema')
    return {name: {"number": number,"generators": generators,
                   "operations": closure(generators)}
            for number,name,generators in data['source_groups']}


def transform(A, R, theta, eta):
    """Neumann real-coefficient action; spatial determinant in 3D equals det(R2) here."""
    d=np.linalg.det(R)
    return (eta if theta else 1)*d*R@A@R.T


def sector_basis(kind):
    if kind=='symmetric':
        return (np.array([[1.,0],[0,0]]), np.array([[0.,0],[0,1.]]),
                np.array([[0.,1],[1.,0]]))
    if kind=='antisymmetric': return (np.array([[0.,1],[-1.,0]]),)
    raise ValueError(kind)


def projector_basis(operations, kind, eta, tolerance=1e-8):
    """Reynolds group-average projector on rank-2 tensor basis."""
    bs=[]
    for B in sector_basis(kind):
        Bp=sum((transform(B,R,t,eta) for R,t in operations),np.zeros((2,2)))/len(operations)
        bs.append(Bp.reshape(-1))
    s=np.linalg.svd(np.column_stack(bs).T,compute_uv=False)
    rank=int(np.sum(s>tolerance))
    return rank, None


def independent_nullity(operations,kind,eta,tolerance=1e-8):
    """Independent direct linear constraint null-space, no averaging."""
    basis=np.stack([B.reshape(-1) for B in sector_basis(kind)],axis=1)
    blocks=[]
    for R,t in operations:
        A=np.column_stack([(transform(B,R,t,eta)-B).reshape(-1)
                           for B in sector_basis(kind)])
        blocks.append(A)
    s=np.linalg.svd(np.vstack(blocks),compute_uv=False)
    rank=int(np.sum(s>tolerance))
    return basis.shape[1]-rank


def allowed(group,sector):
    kind,eta=SECTORS[sector]
    pr,_=projector_basis(group['operations'],kind,eta)
    nu=independent_nullity(group['operations'],kind,eta)
    if pr!=nu: raise AssertionError(f'projector/nullspace disagreement: {sector}: {pr} {nu}')
    return pr


def transverse_symmetric_allowed(group, sector='gN',tol=1e-8):
    """Tests the actual sigma_xy+sigma_yx functional, not generic symmetric rank."""
    kind,eta=SECTORS[sector]
    if kind!='symmetric': raise ValueError('symmetric sector only')
    for B in sector_basis(kind):
        P=sum((transform(B,R,t,eta) for R,t in group['operations']),np.zeros((2,2)))/len(group['operations'])
        if abs(P[0,1]+P[1,0])>tol: return True
    return False


# Literal transcription from the source PDF Table I, page 4 (not corrected).
# All 31 2D magnetic groups in the OmegaN row: 'All'.
SOURCE_PRINTED = {
  'gN': "1 m 2 2' mm2 m'm2' 3 4 6 m'm'2 3m' 4' 4m'm' 4'm'm 6m'm'".split(),
  'gA': "1 m 2 2' mm2 m'm2' 3 4 6 11' m' m1' 21' mm21' m'm'2 31' 41' 4' 61' 6'".split(),
  'omegaA': "1 m 2 2' mm2 m'm2' 3 4 6 m' 3m 4mm 6mm".split(),
}


def audit_table():
    groups=load_groups()
    source={**SOURCE_PRINTED, 'omegaN':list(groups)}
    table=[]
    for name,group in groups.items():
        computed={s:allowed(group,s)>0 for s in SECTORS}
        paper={s:name in source[s] for s in SECTORS}
        table.append({"number":group['number'],"group":name,"order":len(group['operations']),
                      "computed":computed,"paper":paper,
                      "disagreements":[s for s in SECTORS if computed[s]!=paper[s]],
                      "gN_transverse_xy":transverse_symmetric_allowed(group)})
    return table
