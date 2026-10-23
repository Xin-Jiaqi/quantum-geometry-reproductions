#!/bin/bash
# LOCAL-ONLY toolchain environment (git-ignored; do NOT commit).
# Fill in per approved platform. No internal paths belong in committed scripts.
export VASP_BIN=          # e.g. path/to/vasp_ncl (SOC build)
export MPI_RUN=mpirun     # e.g. mpirun (Intel MPI) or mpirun.openmpi
# optional for para VaspWannier build:
# export LD_LIBRARY_PATH=...mkl/intel64:...mpi/lib:...
