"""Same-k PAW overlaps between zero-gate and independent gate-SCF states."""
from pathlib import Path
import os, socket, sys, json, hashlib
import numpy as np
from ase.io import read
if not os.getenv('PBS_JOBID') or not os.getenv('PBS_NODEFILE') or socket.gethostname().startswith(('super','cms-ssh')):
    raise SystemExit('PBS compute node required')
parent=Path('/work/sunq/xjq/mosse-7269-mechanism-20261104')
sys.path[:0]=[str(parent),str(parent/'pawlib')]
from paw import pawpotcar, nonlq
from vaspwfc import vaspwfc
base=Path('/work/sunq/xjq/inse-nmee-20261109/A0/A/seed')
gates=Path('/work/sunq/xjq/inse-nmee-20261109')
root=Path.cwd(); atoms=read(base/'POSCAR')
pots=[pawpotcar(s) for s in (base/'POTCAR').read_text().split('End of Dataset')[:-1]]
w0=vaspwfc(str(base/'WAVECAR'),lsorbit=True)
waves={case:vaspwfc(str(gates/case/'A/seed/WAVECAR'),lsorbit=True) for case in ('Am','Ap')}
for case,w in waves.items():
    a=read(gates/case/'A/seed/POSCAR')
    assert np.max(abs(a.positions-atoms.positions))<1e-9
    assert np.max(abs(a.cell.array-atoms.cell.array))<1e-9
    assert (gates/case/'A/seed/POTCAR').read_bytes()==(base/'POTCAR').read_bytes()
    assert w._nbands==w0._nbands and np.max(abs(w._kvecs-w0._kvecs))<1e-9
def states(w,i):
    return np.stack([w.readBandCoeff(ikpt=i+1,iband=b+1).reshape(2,-1) for b in range(w._nbands)],axis=-1)
out={case:[] for case in waves};normerr=[]
for i,k in enumerate(w0._kvecs):
    g=w0.gvectors(i+1);z=states(w0,i)
    pr=nonlq(atoms,w0._encut,potcar=pots,k=k)
    assert np.array_equal(g,pr.Gvec)
    def overlap(left,right):
        o=sum(left[s].conj().T@right[s] for s in (0,1))
        for ia,typ in enumerate(pr.element_idx):
            op=pr.qproj[typ]*pr.cqfak[typ][:,None]*pr.crexp[None,:,ia]
            q=pots[typ].get_Qij()
            o+=sum((op@left[s]).conj().T@q@(op@right[s]) for s in (0,1))
        return o
    normerr.append(float(abs(overlap(z,z)-np.eye(w0._nbands)).max()))
    for case,w in waves.items():
        assert np.array_equal(w.gvectors(i+1),g)
        f=states(w,i);normerr.append(float(abs(overlap(f,f)-np.eye(w._nbands)).max()))
        out[case].append(overlap(f,z))
    print('CROSS_FIELD_OVERLAP',i+1,len(w0._kvecs),flush=True)
assert max(normerr)<1e-4
np.savez_compressed(root/'overlaps.npz',**{case:np.array(v) for case,v in out.items()})
record={'job':os.environ['PBS_JOBID'],'node':socket.gethostname(),'orientation':'<u_F|u_0>, same k and physical PAW metric',
        'max_norm_error':max(normerr),'zero_source':str(base),'gate_source':str(gates),
        'cache_sha256':{name:hashlib.sha256(p.read_bytes()).hexdigest() for name,p in [
            ('zero',base.parent/'cache.npz'),*[(c,gates/c/'A/cache.npz') for c in waves]]}}
(root/'overlaps.json').write_text(json.dumps(record,indent=2));print('OVERLAPS_COMPLETE',flush=True)
