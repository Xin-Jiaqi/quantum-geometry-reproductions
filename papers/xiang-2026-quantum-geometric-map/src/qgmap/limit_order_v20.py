"""v20: finite regulator -> real-probe-frequency limit of the v19 SCBA Ward kernel.

At fixed NONZERO (q,omega), taking eta->0 removes the 2i*eta*Gamma0
regulator leakage. At (q,omega)=(0,0), Gamma0 diverges ~1/eta and
2*eta*Gamma0 remains finite due to the imaginary disorder self-energy.
This numerical limit is not a full electromagnetic conductivity.
"""
from __future__ import annotations
from dataclasses import replace
import numpy as np
from .warped_disorder_v17 import Lattice
from .finite_qomega_ward_v19 import finite_ward,two_energy_greens


def scan_limit(*,n=61,eta_values=(.08,.02,.005,.002,.001),q=.06,omega=.04):
    if abs(q)<1e-14 and abs(omega)<1e-14:
        raise ValueError('finite probe needs at least one nonzero q or omega')
    entries=[]
    for eta in eta_values:
        p=replace(Lattice(n=n),eta=eta)
        regular=finite_ward(p,q=q,omega=omega)
        origin=finite_ward(p,q=0.,omega=0.)
        x=two_energy_greens(p,qy=q,omega=omega)
        # SCBA scattering exists even as external infinitesimal eta vanishes.
        s=x['sr']
        rates=-np.linalg.eigvalsh((s-s.conj().T)/(2j))
        entries.append(dict(eta=float(eta),fixed_probe_ward_abs=regular['residuals']['wti'],
            fixed_probe_regulator_leak_abs=regular['residuals']['physical_frequency_without_regulator'],
            fixed_probe_gamma_norm=regular['density_vertex_norm'],
            uniform_static_gamma_norm=origin['density_vertex_norm'],
            uniform_static_regulator_leak_abs=origin['residuals']['physical_frequency_without_regulator'],
            uniform_static_ward_abs=origin['residuals']['wti'],
            disorder_scba_broadening_min=float(rates.min()),
            disorder_scba_broadening_max=float(rates.max()),
            fixed_probe_iterations=regular['iterations'],
            zero_probe_iterations=origin['iterations']))
    return dict(protocol='SCBA_regulator_limit_fixed_qomega_versus_uniform_static_origin',
      n=int(n),q=float(q),omega=float(omega),
      field_protocol='onsite',
      values=entries,
      warnings=['finite lattice mesh and noncrossing self-energy; eta->0 only numerically probed',
                'density vertex diverges at q=omega=0, so limits do not commute',
                'no complete electromagnetic Kubo-Streda, nor paper-defined EB PHE'])
