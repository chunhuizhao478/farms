#!/bin/bash
# =============================================================================
# SPE-HFTC 2D section, Gaussian-scattered toughness, eta = 1, full 20-pulse train
# USC CARC Discovery (SLURM). Lines marked <<FILL>> are site/account specific.
# Submit from this folder:  sbatch submit_discovery.sbatch
# Continue after a time-out:  sbatch --export=ALL,RECOVER=1 submit_discovery.sbatch
# =============================================================================
#SBATCH --job-name=spe2d_gauss_eta1
#SBATCH --account=<<FILL: your Discovery project id, e.g. xxx_123>>
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=%x.o%j
#SBATCH --error=%x.e%j
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<<FILL: email>>

module purge
# <<FILL>> the modules your farms build on Discovery was compiled with, e.g.:
# module load gcc/11.3.0 openmpi/4.1.4
FARMS=<<FILL: /path/to/farms-opt on Discovery>>

# --- case definition (everything else comes from ../../elasticity.i) ---------
CASE_DIR="$SLURM_SUBMIT_DIR"
CODE_DIR="$(cd "$CASE_DIR/../.." && pwd)"          # pulsepower/spe_hftc/code_2d
MESH=../mesh/spe_section_h0p3.msh                  # 67 k triangles, l = 0.6 mm (regenerate with ../mesh/make_spe_section_2d.py --h 0.3e-3)
ARGS="mesh_file=$MESH eta=1.0 Np=20 MultiApps/fracture/input_files=fracture_gaussian.i Outputs/file_base=cases/discovery_gaussian_eta1/out"
# Gaussian field knobs live in fracture_gaussian.i (gc_cv = 0.12, gc_min_frac = 0.3, gc_seed = 1); override by editing that
# file or by appending to the sub-app cli_args in elasticity.i. Energy fraction, pulse count, damping: override here, e.g. eta=0.7.

cd "$CODE_DIR" || exit 1
if [[ "${RECOVER:-0}" == "1" ]]; then
    srun --mpi=pmix "$FARMS" -i elasticity.i $ARGS --recover
else
    srun --mpi=pmix "$FARMS" -i elasticity.i $ARGS
fi
# If srun/pmix is not configured for your MPI, use instead:
#   mpirun -np $SLURM_NTASKS "$FARMS" -i elasticity.i $ARGS
