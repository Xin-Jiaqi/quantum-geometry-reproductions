"""Freeze reproducible valley A/B/C coefficients with exact projector jets.

Run from studies/bsf-nmee-connection:
OPENBLAS_NUM_THREADS=1 python scripts/run_stage7e_analytic_coefficients.py
"""
from pathlib import Path
from dataclasses import replace
import json,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from bsf_nmee.sliding_tight_binding import SlidingSpec
from bsf_nmee.projector_jet_response import analytic_valley_ABC


def run():
    records=[]
    for u in (0.,.1,.2,.4):
        spec=replace(SlidingSpec(),shift=(u*.23,u*.41))
        val=analytic_valley_ABC(spec,radius=.75,nr=26,nt=72)
        records.append({"u":u,**val})
        print(f"u {u:.2f}: "+
            " ".join(f"{k}={v['total']:.9g}" for k,v in val["coefficients"].items()),
            flush=True)
    return {
        "scope":"0 K gapped equilibrium/static response, not driven DC",
        "coefficients":"Gamma = A*s^2 + B*c^2 + C*s*c + O(t^4)",
        "integration":"K and Kprime radius 0.75 nm^-1, radial 26, angular 72",
        "rows":records}


if __name__=="__main__":
    result=run()
    p=Path(__file__).resolve().parents[1]/"results/sliding_stage7e_ABC_recomputed.json"
    p.write_text(json.dumps(result,indent=2)+"\n")
    print("Saved",p)
