"""Source-paper model capabilities in one entry, no bulk/device substitution."""
import json
from qgmap.dirac import eq11_kernel,continuum_temperature_kernel
from qgmap.spin_multiband_v5 import four_band_audit
from qgmap.orbital_phe_v4 import orbital_phe_eq56
from qgmap.orbital_mnhe_v5 import orbital_mnhe_2d,orbital_mnhe_3d
from qgmap.orbital_geometry_v54 import band_geometry,periodic_three_band_model,moment_curvature_kernel

def run():
    spin=four_band_audit(N=41)
    phe=orbital_phe_eq56(N=121)
    mnhe2=orbital_mnhe_2d(N=121)
    mnhe3=orbital_mnhe_3d(N=25)
    _,H,V=periodic_three_band_model(N=40)
    moment=moment_curvature_kernel(band_geometry(H,V),mu=-1.1,temperature=.5)
    assert abs(spin['closing_residual'])<1e-6
    assert abs(phe['printed_line_discrepancy'])<1e-6
    assert abs(phe['finite_difference_minus_first'])<1e-6
    assert abs(mnhe2['finite_field_minus_eq48'])<1e-6
    assert abs(mnhe3['finite_field_minus_eq48'])<1e-6
    return dict(spin_PHE={'printed_eq11_kernel':float(eq11_kernel(.125)),
                          'finiteT':float(continuum_temperature_kernel(.02,100.,.125))},
                spin_MNHE=spin,orbital_PHE=phe,orbital_MNHE_2D=mnhe2,orbital_MNHE_3D=mnhe3,
                multiband_moment_sector_only=moment,
                scope='Source-model kernels are callable; general multiband positional shift and physical metallic DC are not inferred',
                source_corrections={'Eq56':'Two lines agree; historical code misread d_ab f as d_ab epsilon*fprime',
                                    'Eq49_to_52':'Factor 2 survives IBP in retained velocity-metric term; source Eq52 prints 1',
                                    'Eq20_to_31':'Relative sign inconsistency remains',
                                    'three_band':'Source Eq49 and Eq52 explicitly drop these; full moment sector is not full MNHE'})

if __name__=='__main__':print(json.dumps(run(),indent=2,allow_nan=False))
