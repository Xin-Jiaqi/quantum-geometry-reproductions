"""Extract W2Se4/Mo2S4/W2S4 six-atom candidates from an original ASE BiDB SQLite file.

Only a SOURCE-PROVENANCE PRE-SCREEN. Does not declare AB/BA or full
structural relaxation. Never changes the frozen Jia et al. reproduction.
Requires ASE when running on a real database.
"""
import argparse
from collections import Counter
import hashlib, json, csv, sqlite3
from pathlib import Path

TARGETS={"WSe2":{"W":2,"Se":4},"MoS2":{"Mo":2,"S":4},"WS2":{"W":2,"S":4}}
METADATA=("bilayer_uid","monolayer_uid","number_of_layers","slide_stability",
          "binding_energy_gs","binding_energy_zscan","distance","gap_pbe",
          "space_group_number","layer_group_number")


def identify(symbols):
    species=Counter(symbols)
    return next((name for name, counts in TARGETS.items()
                 if species==Counter(counts)),None)


def sha256_file(path):
    sha=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            sha.update(chunk)
    return sha.hexdigest()


def ase_sqlite_shape(path):
    con=sqlite3.connect(f"file:{Path(path).resolve()}?mode=ro",uri=True)
    try:
        tables={row[0] for row in con.execute(
            "SELECT name FROM sqlite_master WHERE type='table'")}
        return {"systems","species"}.issubset(tables)
    finally:
        con.close()


def extract(database,out):
    database=Path(database);out=Path(out)
    if not database.is_file() or not ase_sqlite_shape(database):
        raise ValueError("Not a readable original ASE-SQLite db; no geometry accepted")
    from ase.db import connect
    from ase.io import write
    out.mkdir(parents=True,exist_ok=True)
    records=[]
    with connect(str(database)) as db:
        for row in db.select():
            atoms=row.toatoms()
            material=identify(atoms.get_chemical_symbols())
            if not material:continue
            name=f"{material}_bidb_row_{row.id}_UNVERIFIED_GEOMETRY.vasp"
            dest=out/name
            write(str(dest),atoms,format="vasp",direct=True,sort=True)
            record={"row_id":row.id,"formula":material,
                    "filename":name,"export_sha256":sha256_file(dest),
                    "geometry_status":"UNVERIFIED_AB_BA_AND_RELAXATION"}
            for key in METADATA:
                v=row.key_value_pairs.get(key)
                record[key]=v if v is None or isinstance(v,(str,int,float,bool)) else repr(v)
            records.append(record)
    package={"source_sha256":sha256_file(database),
             "status":"GEOMETRY_EXTRACTED_NOT_AUTHENTICATED_OR_RELAXED",
             "rows":records,
             "warning":"No AB/BA stacking assignment or full DFT-relaxation inferred"}
    (out/"candidates.json").write_text(json.dumps(package,indent=2)+"\n")
    with (out/"candidates.csv").open("w",newline="") as f:
        fields=("row_id","formula","filename","export_sha256",
                "geometry_status",*METADATA)
        writer=csv.DictWriter(f,fieldnames=fields,extrasaction="ignore")
        writer.writeheader();writer.writerows(records)
    return records


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("database")
    parser.add_argument("--out",default="extracted_targets")
    p=parser.parse_args()
    result=extract(p.database,p.out)
    print("Extracted",len(result),"unverified bilayer geometry candidates")
