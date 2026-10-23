#!/bin/bash
#SBATCH --job-name=BD36-2255_bilayer
#SBATCH --partition=root.bingxing
#SBATCH --nodes=1
#SBATCH --ntasks=28
#SBATCH --time=02:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
set -eo pipefail
ROOT="${SLURM_SUBMIT_DIR:?}"
export LD_LIBRARY_PATH=/share/home/e6a000183/intel/oneapi/mkl/2022.1.0/lib/intel64:$LD_LIBRARY_PATH
export PATH=/share/home/e6a000183/tools/vasp/VaspWannier/vasp.5.4.4/bin:$PATH
export PATH=/share/home/e6a000183/intel/oneapi/mpi/2021.6.0/bin:$PATH
MPI=mpirun
VASP=vasp_ncl
test -x "$(command -v $MPI)" && test -x "$(command -v $VASP)"
cd "$ROOT"
$MPI -np 28 "$VASP" > vasp.stdout
echo "EXIT_OK"
