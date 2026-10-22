#!/bin/bash
# =============================================================================
# SPE-HFTC 2D section, Gaussian-scattered toughness, eta = 0.3, full 20-pulse train
# USC CARC Discovery (SLURM). Account, modules and farms path filled on 2026-10-21 for user zhaochun.
# Submit from this folder:  sbatch submit_discovery.sbatch
# Continue after a time-out:  sbatch --export=ALL,RECOVER=1 submit_discovery.sbatch
# =============================================================================
#SBATCH --job-name=spe2d_gauss_eta0p3
#SBATCH --account=scec_608
#SBATCH --partition=oneweek      # main had 0 idle nodes and >1000 pending jobs on 2026-10-21; oneweek had 18 idle
#SBATCH --nodes=1
#SBATCH --ntasks=20              # oneweek nodes have 20 cores
#SBATCH --cpus-per-task=1
#SBATCH --mem=58G
#SBATCH --time=24:00:00
#SBATCH --output=%x.o%j
#SBATCH --error=%x.e%j

module purge
module load usc/13.3.0 gcc/13.3.0 openmpi/5.0.5 openblas/0.3.28 python/3.11.9   # the set the login shell loads; farms-opt was built against this MOOSE source tree
FARMS=/home1/zhaochun/projects/farms_cdms/farms-opt

# --- case definition (everything else comes from ../../elasticity.i) ---------
CASE_DIR="$SLURM_SUBMIT_DIR"
CODE_DIR="$(cd "$CASE_DIR/../.." && pwd)"          # pulsepower/spe_hftc/code_2d
MESH=../mesh/spe_section_h0p3.msh                  # 67 k triangles, l = 0.6 mm (regenerate with ../mesh/make_spe_section_2d.py --h 0.3e-3)
ARGS="mesh_file=$MESH eta=0.3 Np=20 MultiApps/fracture/input_files=fracture_gaussian.i Outputs/file_base=cases/discovery_gaussian_eta0p3/out"
# Gaussian field knobs live in fracture_gaussian.i (gc_cv = 0.12, gc_min_frac = 0.3, gc_seed = 1); override by editing that
# file or by appending to the sub-app cli_args in elasticity.i. Energy fraction, pulse count, damping: override here, e.g. eta=0.7.

cd "$CODE_DIR" || exit 1
echo "job $SLURM_JOB_ID on $(hostname), $SLURM_NTASKS tasks, started $(date)"
if [[ "${RECOVER:-0}" == "1" ]]; then
    mpirun -np $SLURM_NTASKS "$FARMS" -i elasticity.i $ARGS --recover
else
    mpirun -np $SLURM_NTASKS "$FARMS" -i elasticity.i $ARGS
fi
echo "finished $(date)"
# Open MPI 5 inside a SLURM allocation; alternative launcher:  srun --mpi=pmix "$FARMS" -i elasticity.i $ARGS
