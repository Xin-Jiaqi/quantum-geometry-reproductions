"""Round16: reproducible generator-to-microscopic 31-group audit (no DFT)."""
import json
import numpy as np
from flux_lab.group_microscopic import (load_groups,local_group_residual,
 verify_catalog_projectors,clean_bz_even_odd,project_frequency_tensor)

def report():
    local={name:local_group_residual(name,[.41,-.32]) for name in load_groups()}
    screens=verify_catalog_projectors()
    bz={}
    for name in ("11'","2'","m'","4'","4mm"):
        even,odd=clean_bz_even_odd(name,n=8)
        bz[name]=dict(mesh=8,max_even=float(np.max(abs(even))),max_odd=float(np.max(abs(odd))),
          even_projection_error=float(np.max(abs(even-project_frequency_tensor(name,even,'even')))),
          odd_projection_error=float(np.max(abs(odd-project_frequency_tensor(name,odd,'odd')))))
    return dict(schema="cui-flux-group-microscopic-round16-v1",date="2026-10-22",
       original_paper="Cui et al. arXiv:2604.09725v3",
       group_count=len(local),local_model_k=[.41,-.32],drive_frequency=[.051,.011],
       catalog_sector_rank_checks=len(screens),rank_checks=screens,
       model_residuals=local,BZ_square_groups=bz,
       scientific_scope="All 31 local nondegenerate synthetic interband Kubo, not full Table I or material KME",
       omitted="non-Abelian bands, collisions, SI, real material, C3/C6 hexagonal full BZ")

if __name__=="__main__":print(json.dumps(report(),ensure_ascii=False,indent=2))
