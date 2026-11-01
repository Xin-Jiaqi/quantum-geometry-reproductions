"""Core IGMC model recipe: C/D, metal Fermi surface, finite omega and SI."""
import json
import numpy as np
from igmc_lab.core import (energy_from_frequency,physical_response_si,dirac_grid,
    dirac_hamiltonian,dirac_response,rashba_analytic,rashba_integrate)
from igmc_lab.fermi_contour import integrate_contour,integrate_conduction_fermi_line
from igmc_lab.green_wavepacket_audit import green_eq44_resolvent,wavepacket_cos_sin,analytical_cosine_harmonics
from igmc_lab.independent_checks import orbital_liouville_kernel


def run():
    omega=energy_from_frequency(1e12,protocol='cyclic_hz')
    rows=[]
    for mu in (0.,1.):
        outputs=[integrate_contour(mu,outer_order=n,inner_order=96,photon_eV=omega) for n in (96,160)]
        a,b=outputs
        fs=integrate_conduction_fermi_line(mu,order=320)
        rows.append(dict(mu_eV=mu,C_xy_A_per_m_T=physical_response_si(b['C_xy']),
            D_xx_A_per_m_T=physical_response_si(omega*b['D_peromega_xx']),
            full_subgap_K_imag_A_per_m_T=physical_response_si(b['K_imag_xx']),
            C_mesh_error=abs(a['C_xy']-b['C_xy']),D_per_omega_mesh_error=abs(a['D_peromega_xx']-b['D_peromega_xx']),
            C_fermi_surface_error=abs(b['C_xy']-fs['C_xy_eq13']),fermi_crossings=b['max_x_crossings']))
    x,y=dirac_grid(101); H,V=dirac_hamiltonian(x,y)
    gf=green_eq44_resolvent(H,V['x'],'x',omega)
    wp=wavepacket_cos_sin(H,V['x'],'x',omega)['K_from_wavepacket']
    liouville=orbital_liouville_kernel(101,omega_eV=omega)
    spectral=dirac_response(101,omega_eV=omega,alpha='x',beta='x')
    full=complex(spectral['Kubo_real'],spectral['Kubo_imag'])
    residuals=dict(GF_wavepacket_per_k=float(np.max(abs(gf-wp))),
                   GF_spectral_integral=float(abs(np.mean(gf)-full)),
                   Liouville_spectral_integral=float(abs(liouville-full)))
    assert max(residuals.values())<3e-12
    assert max(row['C_fermi_surface_error'] for row in rows)<2e-8
    return dict(frequency_Hz=1e12,frequency_protocol='cyclic_hz',photon_eV=omega,g_factor=2,
        Rashba_C_plateau_A_per_m_T=physical_response_si(rashba_analytic(.009)),
        Rashba_numeric_kernel=rashba_integrate(.009,nrad=2400,nphi=48),
        Rashba_analytic_kernel=rashba_analytic(.009),Dirac=rows,independent_residuals=residuals,
        cosine_harmonics=analytical_cosine_harmonics(),
        source_verdict={'SM34':'Printed Q=2Re(r sigma) differs from main Q=Re; do not double physical output',
                        'SM46':'Cosine transform has pi per delta, not printed 2pi; Eq47 restores half harmonic',
                        'Fig2d':'Author raw data unavailable; no inferred origin of amplitude discrepancy',
                        'Eq14':'Conduction/displacement channel parities kept distinct; not a certified whole source convention'},
        scope='Published Rashba/Dirac two-band models and reusable physical-Pauli two-band kernels; no generic scattering/resonant/material claim')

if __name__=='__main__':print(json.dumps(run(),indent=2,allow_nan=False))
