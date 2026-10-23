"""Local-only Round34 source discovery and reproducible privacy-conscious receipt.

READ-ONLY input, NO solver submission, NO network, NO uploads, NO raw file
copy, NO authenticity claim. Limit searches to explicitly selected roots.
Optional mathematical probe of ASCII spn/u/u_dis/eig; never assume the
source-provenance/gauge link without verifying original win/wout/DFT.
"""
from __future__ import annotations
import argparse,hashlib,json,os,sys
from pathlib import Path
from datetime import datetime,timezone

SKIP_DIR={".git",".venv","venv","node_modules",".cache","__pycache__",
          ".tox",".next","site-packages","tmp","scratch"}
SUFFIXES=[
 ("spn_fmt",".spn.fmt"),("u_dis","_u_dis.mat"),("u","_u.mat"),
 ("v","_v.mat"),("hr","_hr.dat"),("r","_r.dat"),
 ("eig",".eig"),("win",".win"),("wout",".wout"),
 ("spn_original",".spn"),("nnkp",".nnkp")]
REQUIRED=["spn_fmt","u","eig","hr","r","win","wout"]

def digest(path,max_bytes):
    if path.stat().st_size>max_bytes:return None
    h=hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b""):
            h.update(chunk)
    return h.hexdigest()

def collect(roots,max_files=60000,max_depth=7,max_size_mib=512):
    groups={};nfiles=0;cutoff=False
    for root in roots:
        p=Path(root).expanduser().resolve()
        if not p.is_dir() or p==Path("/") or p==Path.home():
            raise ValueError("provide existing bounded research-data directories, not / or entire home")
        for directory,dirs,files in os.walk(p,followlinks=False):
            current=Path(directory)
            depth=len(current.relative_to(p).parts)
            dirs[:]=[d for d in dirs if d not in SKIP_DIR and depth<max_depth and
                     not (current/d).is_symlink()]
            for name in sorted(files):
                nfiles+=1
                if nfiles>max_files:
                    cutoff=True;break
                found=next(((role,name[:-len(suff)]) for role,suff in SUFFIXES
                            if name.endswith(suff)),None)
                if found is None:continue
                role,seed=found
                path=current/name
                if not path.is_file() or path.is_symlink():continue
                if role=="spn_original":
                    # Actual formatted seed.spn may be ASCII; some sources are Fortran binary.
                    try:
                        with path.open("rb") as fh:header=fh.read(512).decode("ascii")
                        lines=header.splitlines()
                        maybe=(len(lines)>=2 and len(lines[1].split())==2
                               and all(int(x)>0 for x in lines[1].split()))
                    except (UnicodeError,ValueError,IndexError):
                        maybe=False
                    role="spn_fmt" if maybe else "spn_binary"
                key=hashlib.sha256(str(current/seed).encode()).hexdigest()[:16]
                group=groups.setdefault(key,dict(seed_hash=key,present={},
                            source_directory_disclosed=False))
                size=path.stat().st_size
                group["present"][role]=dict(size_bytes=size,
                    sha256=digest(path,int(max_size_mib*1024**2)),
                    hash_skipped_oversize=size>int(max_size_mib*1024**2))
            if cutoff:break
        if cutoff:break
    summary=[]
    for key,group in sorted(groups.items()):
        found=set(group["present"])
        missing=sorted(set(REQUIRED)-found)
        if "spn_fmt" not in found and "spn_binary" in found:
            missing.append("FORMATTED_SPN_CONVERSION_REQUIRED_OFFLINE")
        if "u_dis" not in found:
            missing.append("U_DIS_REQUIRED_WHEN_DISENTANGLED_VALIDATE_NBANDS")
        group["missing_or_unverified"]=sorted(set(missing))
        group["plausible_ascii_chain"]=not(set(REQUIRED)-found)
        group["gauge_and_material_provenance_authenticated"]=False
        group["can_claim_real_kme"]=False
        group["physical_spin_source_authenticated"]=False
        summary.append(group)
    return dict(schema="cui-round34-local-source-collection-v1",
       paper_status="PARTIAL",real_material_KME_verified=False,
       original_files_uploaded=False,computations_submitted=False,
       scanned_file_count=nfiles,search_cutoff_reached=cutoff,
       found_seed_groups=len(summary),
       groups=summary)

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root",action="append",default=[],
        help="explicit research source folder; may repeat; never entire home")
    parser.add_argument("--out",default="local_evidence/round34/receipt.json")
    parser.add_argument("--max-files",type=int,default=60000)
    parser.add_argument("--max-depth",type=int,default=7)
    parser.add_argument("--max-hash-mib",type=int,default=512)
    args=parser.parse_args(argv)
    roots=args.root or [x for x in os.environ.get("CUI_SOURCE_ROOTS","").split(os.pathsep) if x]
    if not roots:
        parser.error("specify --root or CUI_SOURCE_ROOTS; never infer a private data path")
    if not (1<=args.max_files<=100000 and 1<=args.max_depth<=12 and
            1<=args.max_hash_mib<=4096):
        parser.error("scan bounds outside safe range")
    report=collect(roots,args.max_files,args.max_depth,args.max_hash_mib)
    dest=Path(args.out).expanduser()
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(dict(receipt=str(dest),group_count=report["found_seed_groups"],
       scan_cutoff=report["search_cutoff_reached"],
       complete_ascii_candidate_count=sum(g["plausible_ascii_chain"] for g in report["groups"]),
       all_authentication_unverified=True,paper="PARTIAL"),ensure_ascii=False))

if __name__=="__main__":main()
