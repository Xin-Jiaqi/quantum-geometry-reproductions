"""Extract genuine BiDB z-scan geometry records without installing ASE.

Original SQLite ASE arrays are little-endian: Z=int32, xyz/cell=float64.
Preserve source ID/hash and NEVER label records fully relaxed or DFT-NMEE.
Run: python material-screening/extract_bidb_without_ase.py bidb.db --out /tmp/tmd
"""
from collections import Counter
from pathlib import Path
import argparse, hashlib, importlib.util, json, sqlite3
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
checker=ROOT/"material-screening"/"ideal_stack_symmetry.py"
modspec=importlib.util.spec_from_file_location("stage8b_point_checker",checker)
checker_module=importlib.util.module_from_spec(modspec)
modspec.loader.exec_module(checker_module)
CHEM={74:"W",34:"Se",42:"Mo",16:"S"}
TARGETS={(74,34):"WSe2",(42,16):"MoS2",(74,16):"WS2"}

def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for part in iter(lambda:f.read(1<<20),b""):h.update(part)
    return h.hexdigest()

def get_record(con,row_id):
    record=con.execute("SELECT natoms,numbers,positions,cell,key_value_pairs FROM systems WHERE id=?",
                       (row_id,)).fetchone()
    if record is None:raise KeyError(row_id)
    n,nb,pb,cb,meta=record
    numbers=np.frombuffer(nb,dtype="<i4").copy()
    cart=np.frombuffer(pb,dtype="<f8").reshape(-1,3).copy()
    cell=np.frombuffer(cb,dtype="<f8").reshape(3,3).copy()
    if len(numbers)!=n or len(cart)!=n:
        raise ValueError("Corrupt ASE atomic arrays")
    return numbers,cart,cell,json.loads(meta)

def assess(numbers,cart,cell):
    frac=cart@np.linalg.inv(cell)
    site=frac.copy()
    site[:,2]=cart[:,2]-np.mean(cart[:,2])
    labels=[CHEM[int(z)] for z in numbers]
    ops=checker_module.operations(site,labels,cell[:2,:2])
    projectors=[]
    for i in range(2):
      for j in range(2):
        M=np.zeros((2,2));M[i,j]=1
        avg=sum(f*R@M@R.T for n,mir,f,R in ops)/len(ops)
        projectors.append(avg.reshape(4))
    rank=int(np.linalg.matrix_rank(np.stack(projectors),tol=1e-8))
    return frac,ops,rank

def export(database,out):
    database=Path(database).resolve();out=Path(out);out.mkdir(parents=True,exist_ok=True)
    con=sqlite3.connect("file:"+str(database)+"?mode=ro",uri=True)
    if con.execute("PRAGMA quick_check").fetchone()[0]!="ok":
        raise ValueError("BiDB SQLite integrity check failed")
    records=[]
    for (zm,zx),material in TARGETS.items():
        matches=con.execute(
            "SELECT s.id FROM systems s JOIN species m ON m.id=s.id AND m.Z=? AND m.n=2 "
            "JOIN species x ON x.id=s.id AND x.Z=? AND x.n=4 WHERE s.natoms=6 ORDER BY s.id",
            (zm,zx)).fetchall()
        for (row_id,) in matches:
            numbers,xyz,cell,meta=get_record(con,row_id)
            if meta.get("number_of_layers")!=2 or Counter(numbers)!=Counter({zm:2,zx:4}):
                continue
            frac,ops,rank=assess(numbers,xyz,cell)
            if len(ops)==6:group="C3v"
            elif len(ops)==12 and any(n==0 and not mir and f==-1 for n,mir,f,R in ops):
                group="D3h"
            elif len(ops)==12 and any(n==3 and not mir and f==-1 for n,mir,f,R in ops):
                group="D3d"
            else:group="OTHER_ORDER_"+str(len(ops))
            name=f"{material}_BiDB_row{row_id}_{group}_zscan_notfullyrelaxed.vasp"
            unique=sorted(set(numbers),key=lambda z:list(numbers).index(z))
            indices=[i for z in unique for i in range(len(numbers)) if numbers[i]==z]
            with (out/name).open("w") as fp:
                fp.write(f"ASE BiDB row {row_id} {meta.get('bilayer_uid')}; zscan NOT full relaxation\n1.0\n")
                for x in cell:fp.write(" ".join(f"{a:.12f}" for a in x)+"\n")
                fp.write(" ".join(CHEM[int(z)] for z in unique)+"\n")
                fp.write(" ".join(str(int(sum(numbers==z))) for z in unique)+"\nDirect\n")
                for i in indices:fp.write(" ".join(f"{a:.12f}" for a in frac[i]%1)+"\n")
            records.append({"id":row_id,"material":material,"bilayer_uid":meta.get("bilayer_uid"),
                "source_slide_stability":meta.get("slide_stability"),
                "source_layer_group":meta.get("layer_group"),
                "source_gap_pbe_eV":meta.get("gap_pbe"),
                "source_interface_gap_A":meta.get("distance"),
                "independent_group":group,"inplane_tensor_dimension":rank,
                "file":name,"file_sha256":digest(out/name),
                "status":"BIDB_ZSCAN_FROZEN_MONOLAYER_NOT_FULLY_RELAXED"})
    outdata={"source_sha256":digest(database),"row_count":len(records),
             "provenance_warning":"Source hash does not independently establish author release hash.",
             "not_calculated":"fully relaxed bilayer, full-path SOC gap, real-material NMEE, driven DC",
             "records":records}
    (out/"source_record_audit.json").write_text(json.dumps(outdata,indent=2)+"\n")
    con.close()
    return outdata

if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("database")
    parser.add_argument("--out",required=True);args=parser.parse_args()
    result=export(args.database,args.out)
    print("Exported",result["row_count"],"real-source z-scan candidates, not fully relaxed")
