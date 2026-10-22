"""Regenerate concise independent finite-q periodic-spin magnetic response evidence."""
import json,argparse
from pathlib import Path
from dataclasses import replace
import numpy as np
from qgmap.bulk_finite_q_spin_v72 import (Spec,spectral,liouville_independent,
    analytic_zeeman_mixed,field_difference,incorrect_identity_vertex,incorrect_local_velocity_vertex)

def pair(z):return {'real':float(z.real),'imag':float(z.imag)}

def generate():
    s=Spec()
    r=spectral(s);m=analytic_zeeman_mixed(s);fd=field_difference(s)
    g=spectral(replace(s,orbital_gauge=.7))
    bad=incorrect_identity_vertex(replace(s,orbital_gauge=.7))
    wrong=incorrect_local_velocity_vertex(s)
    indep=liouville_independent(s)
    return {'version':'v72','source_status':'PARTIAL',
      'scope':'finite-q periodic Bloch spin Zeeman scalar-potential susceptibility; NOT source-selected tau pole or orbital magnetotransport',
      'inputs':{k:getattr(s,k) for k in s.__dataclass_fields__},
      'rho':pair(r['density']), 'jy':pair(r['longitudinal']),
      'jx_per_E':pair(r['electric_normalized_transverse']),
      'dB_jx_per_E':pair(m['dB_E_transverse']),
      'dB_jx_per_E_five_point':pair(fd['electric_normalized_transverse']),
      'independent_dense_jx_per_E':pair(indep['electric_normalized_transverse']),
      'density_ward':pair(r['ward_residual']),
      'magnetic_ward':pair(m['dB_ward_residual']),
      'gauged_jx_per_E':pair(g['electric_normalized_transverse']),
      'wrong_density_vertex_ward_defect':pair(bad['ward_defect']),
      'wrong_endpoint_velocity_ward_defect':pair(wrong['ward_defect'])}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',default='results/frozen/v72_bulk_finite_q_spin.json')
    args=p.parse_args();f=Path(args.output);f.parent.mkdir(parents=True,exist_ok=True)
    f.write_text(json.dumps(generate(),sort_keys=True,indent=2)+'\n',encoding='utf-8')
    print('v72 finite-q evidence',f)
