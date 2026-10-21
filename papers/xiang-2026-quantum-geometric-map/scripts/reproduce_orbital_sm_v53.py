#!/usr/bin/env python3
"""Original PRB SM pp11-12 Eq48-49-52-56 reproducible numerical evidence."""
import argparse, hashlib, json
from pathlib import Path
from qgmap.orbital_sm_v53 import full_audit,SM_PDF_SHA256,MAIN_PDF_SHA256

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sm-pdf',type=Path,help='Optional exact-byte source SHA256 check')
    parser.add_argument('--main-pdf',type=Path)
    parser.add_argument('--output',type=Path)
    a=parser.parse_args()
    for f,sha in ((a.sm_pdf,SM_PDF_SHA256),(a.main_pdf,MAIN_PDF_SHA256)):
        if f is not None and hashlib.sha256(f.read_bytes()).hexdigest()!=sha:
            parser.error('Original source PDF SHA mismatch: '+str(f))
    data=full_audit()
    data['source_pdf_bytes_checked']={'sm':a.sm_pdf is not None,'main':a.main_pdf is not None}
    out=json.dumps(data,indent=2,allow_nan=False)+'\n'
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(out,encoding='utf8')
    else:print(out,end='')
if __name__=='__main__':main()
