"""Independent core reproducer for Xiang et al., PRB 113 L201406 (2026)."""
from .dirac import Parameters, hamiltonian, energies, eq11_kernel, local_zbc_analytic, local_zbc_eigensolver
from .rashba import zero_temperature, direct_fermi_ring_sum
from .liouville import two_color_response, equilibrium_population_derivative
