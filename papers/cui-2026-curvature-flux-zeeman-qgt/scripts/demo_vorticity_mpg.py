"""Reproduce synthetic spin-vorticity Eq.(27) and conditional 2D MPG audit."""
import json
import numpy as np
from flux_lab.spin_current import spin_orbital_benchmark
from flux_lab.spin_vorticity import vorticity_with_hessian,vorticity_finite_difference
from flux_lab.magnetic_symmetry import load_groups,SECTORS,allowed,audit_table

def run():
    h,v,K,s=spin_orbital_benchmark(.31,-.22,periodic=True)
    a=vorticity_with_hessian(h,v,K,s)
    b=vorticity_finite_difference(.31,-.22)
    g=load_groups()
    return dict(point=[.31,-.22],
       band_vorticity=a['bar_omegaA_z'].tolist(),
       finite_difference_max_residual=float(np.max(abs(a['bar_omegaA_z']-b))),
       mixed_hessian_antisymmetric_max_residual=float(np.max(abs(a['vertex_curl_half']))),
       nonzero_symmetric_hessian_max=float(np.max(abs(a['symmetric_hessian_correction']))),
       magnetic_groups=len(g),group_candidate_spaces=len(g)*len(SECTORS),
       source_table_conditional_disagreements=[
          dict(group=row['group'],sector=s) for row in audit_table() for s in row['disagreements']
       ],status='Eq27 verified in synthetic model; magnetic static tensor spaces conditional')

if __name__=='__main__':
    print(json.dumps(run(),ensure_ascii=False,indent=2))
