"""Recompute Stage7j pure orbital-gauge stress tests.

From studies/bsf-nmee-connection:
 OPENBLAS_NUM_THREADS=1 python scripts/run_stage7j_wannier_gauge.py

The independent Stage7j reference JSON is immutable; output is separate.
No actual DFT or Wannierization takes place here.
"""
from pathlib import Path
from dataclasses import replace
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from bsf_nmee.sliding_tight_binding import SlidingBilayer,SlidingSpec
from bsf_nmee.wannier_gauge_links import gauge_stress_case


def run():
    rows=[]
    for u,t in ((.1,.005),(.4,.005),(.4,.020)):
        spec=replace(SlidingSpec(),shift=(u*.23,u*.41),
                     tunneling_eV=t)
        k=np.array([1/3,2/3])@SlidingBilayer(spec).model().reciprocal
        k=k+np.array([.05,.02])
        row=gauge_stress_case(spec,k)
        rows.append(row)
        print(f"u={u} t={t} original={row['original']:.9g} "
              f"corrected={row['corrected']:.9g} "
              f"wrong={row['naive']:.9g} "
              f"physical_err={row['physical_gauge_error']:.2e}",flush=True)
    return {"scope":"8-band synthetic, 0-K insulating STATIC mixed Berry curvature",
            "change":"pure k/lambda dependent orbital representation, no physical perturbation",
            "locality":"first K valley center + (0.05,0.02) inverse nm ONLY",
            "rows":rows,
            "DFT":"NOT TESTED: this is a synthetic gauge gate, not real Wannier data"}


if __name__=="__main__":
    record=run()
    path=ROOT/"results/sliding_stage7j_gauge_recomputed.json"
    path.write_text(json.dumps(record,indent=2)+"\n")
    print("Saved",path)
