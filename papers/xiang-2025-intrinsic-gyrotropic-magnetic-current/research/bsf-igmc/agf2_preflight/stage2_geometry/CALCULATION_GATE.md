# Stage 2: physical eligibility and acceptance gate

**G0 (must pass before any real DFT)**: Valid VASP license/executable; hash-locked Ag/F PAW metadata and private checksum; explicit approval of core-hour budget and actual scheduler; documented pseudopotential order, PBE, Ag-4d Dudarev U_eff=4 eV, DFT-D3 variant, ENCUT=500 eV, electronic convergence 1e-6 eV, ionic forces < 0.01 eV/Angstrom, 20 Angstrom vacuum. kmesh, smearing, MAGMOM, force constraints and dipole flags are *operator choices*, not reported paper parameters. Never commit POTCAR.

**G1 (initial minimal pilot)**: First optimize Type-I B1 from 2.6/3.0/3.4 Angstrom vertical interface-gap guesses; 0% strain, fixed in-plane lattice, collinear spin-polarized non-SOC. Record all attempted runs, raw OUTCAR/vasprun hashes or secure storage paths, magnetic moments, electronic convergence, forces, final total energies and relaxed structures. If spin order collapses, SCF fails or contacts become unphysical, stop and label failures.

**G2 (symmetry and sliding stability)**: A converged B1 minimum may be mapped to B2 using unitary m_z and AXIAL sign-flipped in-plane moments. Verify unique same-element site mapping, energies and forces with identical settings. Relax in-plane slide and inspect an energy cut along that coordinate. **Only a verified local energy double-well** plus polarization can justify sliding-ferroelectric language; matching mirror partners alone cannot.

**G3 (not in this pilot)**: SOC magnetic anisotropy, Type-I vs Type-II competition (published ~0.26 meV/f.u. requires high precision), finite-field Berry polarization, CI-NEB, Wannier operators and IGMC must receive separate scientific and resource clearance.

For result submission: preserve source SHA256 and complete INCAR/KPOINTS/POSCAR, VASP/PAW versions, job scripts, calculations with convergence evidence, relaxed geometry, symmetry, magnetic-order comparison, failed runs and an explicit verdict. `data_origin=NOT_RUN` must not change without raw computational evidence.
