"""Stage 8e: source-backed BiDB POSCAR geometry tests (NO ASE or DFT).

These committed fixtures are extracted from the USER-UPLOADED BiDB format ASE
SQLite file; row IDs and geometries correspond to its original 2-layer
z-scan records, NOT authenticated fully ionically relaxed bilayers.

Reuses the Stage8b *independent* 6-atom symmetry matcher, without importing
the Stage8e source auditor. All calculations use NumPy only.
"""
import importlib.util
import unittest
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
CHECKER=ROOT/"material-screening"/"ideal_stack_symmetry.py"
spec=importlib.util.spec_from_file_location("stage8b_original_checker",CHECKER)
geom=importlib.util.module_from_spec(spec)
spec.loader.exec_module(geom)
DIR=ROOT/"material-screening"/"source-structures"
FIXTURES={
 "WSe2_polar":("WSe2_BiDB_row2685_3R_polar_zscan.vasp",6,1,True),
 "WSe2_AA":("WSe2_BiDB_row2688_AA_reference_zscan.vasp",12,0,False),
 "WSe2_2H":("WSe2_BiDB_row2684_2H_inversion_zscan.vasp",12,0,False),
 "MoS2_polar":("MoS2_BiDB_row1930_3R_polar_zscan.vasp",6,1,True)
}

def load_vasp(path):
    lines=path.read_text().splitlines()
    assert "zscan" in lines[0] and "NOT full relaxation" in lines[0]
    cell=np.array([list(map(float,line.split())) for line in lines[2:5]])
    labels=[name for name,n in zip(lines[5].split(),map(int,lines[6].split())) for _ in range(n)]
    f=np.array([list(map(float,line.split())) for line in lines[8:]])
    assert len(f)==6 and len(labels)==6
    f[:,2]=f[:,2]*cell[2,2]
    f[:,2]-=f[:,2].mean()
    return f,labels,cell

def permitted_tensor_dimension(ops):
    basis=[]
    for i in range(2):
      for j in range(2):
        source=np.zeros((2,2));source[i,j]=1
        avg=np.zeros((2,2))
        for n,refl,zflip,R in ops:
            avg+=zflip*(R@source@R.T)
        basis.append((avg/len(ops)).reshape(-1))
    return np.linalg.matrix_rank(np.stack(basis),tol=1e-8)

class SourceBackedBilayerSymmetryTests(unittest.TestCase):
 def test_full_atom_point_groups_tensor_permission(self):
    for name,(filename,n,dim,allowed) in FIXTURES.items():
        f,labels,cell=load_vasp(DIR/filename)
        ops=geom.operations(f,labels,cell[:2,:2])
        self.assertEqual(len(ops),n,name)
        self.assertEqual(permitted_tensor_dimension(ops),dim,name)
        self.assertEqual(geom.allows_xzx(ops),allowed,name)

 def test_real_wse2_parallel_3r_metal_registry(self):
    f,labels,cell=load_vasp(DIR/FIXTURES["WSe2_polar"][0])
    metal=np.flatnonzero(np.array(labels)=="W")
    shift=(f[metal[1],:2]-f[metal[0],:2])%1.
    np.testing.assert_allclose(shift,[1/3,2/3],atol=1e-7)
    self.assertAlmostEqual(np.linalg.norm(cell[0]),3.319006624468,places=9)

 def test_wse2_aa_reference_is_not_polar_and_2h_inverted(self):
    f,labels,cell=load_vasp(DIR/FIXTURES["WSe2_AA"][0])
    ops=geom.operations(f,labels,cell[:2,:2])
    self.assertTrue(any(n==0 and not mirror and sz==-1 for n,mirror,sz,R in ops))
    f,labels,cell=load_vasp(DIR/FIXTURES["WSe2_2H"][0])
    ops=geom.operations(f,labels,cell[:2,:2])
    self.assertTrue(any(n==3 and not mirror and sz==-1 for n,mirror,sz,R in ops))

if __name__=="__main__":
    unittest.main()
