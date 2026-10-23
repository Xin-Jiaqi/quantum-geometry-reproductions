#!/bin/bash
# Slurm draft (para). NOT_RUN - explicit user approval required.
# Toolchain locations come from the LOCAL env file (run_env.sh, NOT committed).
#SBATCH --job-name=BD36-2255_bilayer
#SBATCH --partition=root.bingxing
#SBATCH --nodes=1
#SBATCH --ntasks=28
#SBATCH --time=2:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
set -eo pipefail
ROOT="${SLURM_SUBMIT_DIR:?}"
cd "$ROOT"
[ -f ./run_env.sh ] || { echo "missing run_env.sh (local, uncommitted)"; exit 2; }
source ./run_env.sh
: "${VASP_BIN:?}"
$MPI_RUN -np 28 "$VASP_BIN" > vasp.stdout
echo "EXIT_OK"
