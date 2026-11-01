"""Compact geometry -> source response -> independent dynamics recipe."""
import json
import numpy as np
from flux_lab.core import (model,qgt_eigen,qgt_projector,flux_circle_from_qgt,
                           winding_ellipse,berry_phase_wilson_circle)
from flux_lab.causal_response import (interband_kubo,interband_liouville,
    source_kme_expansion,source_four_sector_geometric,static_subtracted_gyrotropic,reciprocal_residual)
from flux_lab.antiunitary_response import operator_residual,kubo_symmetry_residual
from demo_causal import synthetic_spin_orbital,bz_example


def run():
    d,dk=model(.37,-.52,.3)
    geometry=float(np.max(abs(qgt_eigen(d,dk)['T']-qgt_projector(d,dk)['T'])))
    H,V,S=synthetic_spin_orbital(.37,-.29); opts=dict(mu=2.3,temperature=.025)
    spectral=interband_kubo(H,V,S,z=.035+.003j,**opts)
    orbital=interband_liouville(H,V,S,z=.035+.003j,**opts)
    dynamic=max(float(np.max(abs(a-b))) for a,b in zip(spectral,orbital))
    kme=[]; ghe=[]
    for w in (.04,.02,.01):
        a=source_kme_expansion(H,V,S,omega_source=w,**opts)
        exact=interband_liouville(H,V,S,z=-w,kappa=.5,**opts)[1]-a['static']
        kme.append(dict(omega_source=w,appendix126_error=float(np.max(abs(exact-a['appendix126_dynamic']))),
                        main19_error=float(np.max(abs(exact-a['main19_printed_dynamic'])))))
        b=source_four_sector_geometric(H,V,S,omega_source=w,**opts)
        c=static_subtracted_gyrotropic(H,V,S,omega_drive=-w,**opts)
        ghe.append(float(np.max(abs(b['sigma_source_omega']-c))))
    grids=[bz_example(n,temperature=.2) for n in (32,64,96)]
    return dict(geometry_eigen_projector_error=geometry,
        massless_flux=flux_circle_from_qgt(1.,points=1024),
        winding=winding_ellipse(1.,1.,points=1024),berry_phase=berry_phase_wilson_circle(1.,points=1024),
        spectral_orbital_error=dynamic,source_GHE_expansion_errors=ghe,source_KME=kme,
        reciprocal_error=reciprocal_residual(H,V,S,omega_drive=.035,eta=.003,**opts)['max_residual'],
        mirror_prime_operator=operator_residual('mxT',np.array([.37,-.52])),
        mirror_prime_response={k:v for k,v in kubo_symmetry_residual('mxT',np.array([.37,-.52])).items() if k in ('chi','alpha')},
        four_band_BZ=grids,
        source_verdict={'Eq16_vs_112':'Agree after static subtraction and omega_source=-omega_drive',
                        'Eq19_vs_126':'Linear terms disagree in the SAME source frequency; independent dynamics supports Eq126',
                        'TableI':'Invalid mirror-prime and universal C3 radial counterexamples rejected; whole table not certified'},
        scope='Clean interband model core, physical Pauli spin; source KME is sigma/2, no SI or full metallic collision claim')

if __name__=='__main__':print(json.dumps(run(),indent=2,allow_nan=False))
