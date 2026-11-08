"""Use the user's point-operation catalog for spatial NMEE selection rules.

No optical j/k symmetrization: electric and magnetic fields are distinct inputs.
Only spatial constraints, not a microscopic model or a thermodynamic identity.
"""
from pathlib import Path
import argparse,json,sys,subprocess
import numpy as np
parser=argparse.ArgumentParser();parser.add_argument('toolkit',type=Path);args=parser.parse_args()
sys.path.insert(0,str(args.toolkit/'src'))
from group_theory_operations import point_group_operations,equivariant_map_basis
records={}
for pg in ('C3v','D3h','D3d','Ci','C2v'):
    ops=point_group_operations(pg);R=[np.array(o.matrix()) for o in ops]
    axial=[np.linalg.det(r)*r for r in R]
    inputs=[np.kron(r,a).tolist() for r,a in zip(R,axial)]
    for layered in (False,True):
        # These layer groups preserve the z axis as an unoriented line.
        assert all(np.max(abs(r[2,:2]))<1e-10 for r in R)
        outputs=[np.kron(np.eye(2) if r[2,2]>0 else [[0,1],[1,0]],a).tolist()
                 if layered else a.tolist() for r,a in zip(R,axial)]
        basis=np.array(equivariant_map_basis(outputs,inputs))
        selected=[]
        for b in basis:
            ten=b.reshape(2 if layered else 1,3,3,3)
            c=ten[:,:2,:2,2]
            if abs(c).max()>1e-9:selected.append(c.tolist())
        records[pg+('_layered' if layered else '_whole')]={'dimension_full':len(basis),'inplane_E_Bz_nonzero_basis':selected}
        for b in basis:
            for a,inp in zip(outputs,inputs):
                assert np.max(abs(np.array(a)@b-b@np.array(inp)))<1e-8
        if not layered and pg in ('D3h','D3d','Ci'):assert not selected
        if pg=='C3v':
            for c in selected:
                c=np.array(c);assert np.max(abs(c[:,0,0]-c[:,1,1]))<1e-9
                assert np.max(abs(c[:,0,1]))<1e-9 and np.max(abs(c[:,1,0]))<1e-9
records['mirror_domain_relation']={'operation':'m001; layers exchanged',
    'formula':'C_lower^B(Fz) = -C_upper^A(-Fz); C_total^B(Fz) = -C_total^A(-Fz)',
    'requires':'mirror-related complete environments; symmetric free bilayer or simultaneous reversal of environmental asymmetry',
    'does_not_guarantee':'existence of a compensation field or sign of an individual layer'}
records['provenance']={'toolkit_commit':subprocess.check_output(['git','-C',str(args.toolkit),'rev-parse','HEAD'],text=True).strip(),
    'equation':'axial output spin = tensor * (polar E tensor-product axial B)',
    'note':'spatial rank-three polar transformation; no E/B exchange symmetry imposed'}
(Path(__file__).parent/'symmetry.json').write_text(json.dumps(records,indent=2))
print(json.dumps(records,indent=2))
