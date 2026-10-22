"""Round21 complement: exact occupancy sector decomposition of clean Kubo.

Do not alter the accepted Round20 bz_hotspots.py interface.
"""
import json
from flux_lab.occupation_partition import pair_resolved,integrate_fermi_sectors
from flux_lab.hexagonal_bz import hexagonal_symmetrized_model

def run():
    rows=[
        integrate_fermi_sectors("6",domain="hexagon",order=8,component=(0,0)),
        integrate_fermi_sectors("6m'm'",domain="hexagon",order=8,component=(0,1)),
        integrate_fermi_sectors("6",domain="primitive",method="gauss",
                                order=12,component=(0,0))
    ]
    H,V,S=hexagonal_symmetrized_model([.5145,-1.0456],"6")
    scan=[]
    for T in (.08,.17,.25):
        d=pair_resolved(H,V,S,temperature=T)
        z=d["components"]["same_side"][0,0]
        scan.append({"temperature":T,"same_side_xx":[z.real,z.imag],
                     "min_intercluster_gap":d["min_intercluster_gap"]})
    return {"schema":"cui-flux-round21-thermal-pair-sector-v1",
            "status":"pair decomposition VERIFIED on models; BZ convergence UNRESOLVED",
            "target":"supplement bz_hotspots with across-mu versus same-side thermal occupations",
            "rows":rows,"fixed_k_temperature_scan":scan,
            "no_double_count":"resonance and small-gap subsets overlap and are not independent response channels"}
if __name__=="__main__":print(json.dumps(run(),indent=2,ensure_ascii=False))
