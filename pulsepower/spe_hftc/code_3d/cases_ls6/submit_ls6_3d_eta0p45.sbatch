#!/bin/bash
# =============================================================================
# SPE-HFTC 3D cube on TACC Lonestar6, geometry of 2026-11-17: 7.94 mm (5/16 in.) borehole and electrode gap centred at
# mid-height (manuscript Table 2, Fig. 1, Fig. 2), mesh spe_cube_d7p94_h0p3.msh. 20-pulse train from t = 0, homogeneous
# Gc, eta = 0.45, absorbing side faces, dt_quiet = 5e-7 through EXTRA. (The runs up to 2026-11-16, case
# cube_homog_eta0p45_absorb, used an 8.3 mm borehole and the gap at Lz/3.)
# 18 nodes x 128 cores = 2304 ranks. The mesh split is built by the job if missing. The split, outputs, run.log and
# checkpoints (about 54 GB each) go to $SCRATCH/spe3d, never to the checkout (project rule: results on scratch only).
# Checkpoints every 100 steps, 2 kept. The job exits with the farms exit code. run_chain_ls6.sh drives the segments
# (stall guard, next segment queued behind the running one).
# Submit from this folder:  sbatch --export=ALL,EXTRA="dt_quiet=5e-7" submit_ls6_3d_eta0p45.sbatch
# Continue:                 sbatch --export=ALL,RECOVER=1,EXTRA="dt_quiet=5e-7" submit_ls6_3d_eta0p45.sbatch
# =============================================================================
#SBATCH -J spe3d_eta0p45
#SBATCH -p normal
#SBATCH -N 18
#SBATCH -n 2304
#SBATCH -t 48:00:00
#SBATCH -A EAR24034
#SBATCH -o spe3d_eta0p45.o%j
#SBATCH -e spe3d_eta0p45.e%j

source $WORK/spe3d/farms_cdms/pulsepower/spe_hftc/code_3d/cases_ls6/env_ls6.sh
source $FARMS_DIR/pulsepower/spe_hftc/code_3d/cases_ls6/guard.sh
FARMS=$FARMS_DIR/farms-opt
CODE=$FARMS_DIR/pulsepower/spe_hftc/code_3d
cd $CODE || exit 1
NR=$SLURM_NTASKS
MESH=../mesh/spe_cube_d7p94_h0p3.msh
SPLIT=$SCRATCH/spe3d/spe_cube_d7p94_h0p3_split${NR}.cpr
COMMON="mesh_file=$MESH eta=0.45 zeta_pulse=1e-8 Np=20 ${EXTRA:-}"
CASE=cube_d7p94_mid_eta0p45_absorb
REC=""; [[ "${RECOVER:-0}" == "1" ]] && REC="--recover"
OUT=$SCRATCH/spe3d/$CASE
mkdir -p $OUT || exit 2
[ -e $MESH ] || { echo "ERROR: $MESH missing (copy it into the checkout, it is not in git)"; exit 2; }
if [ ! -d $SPLIT ]; then
  echo "[$(date)] splitting $MESH for $NR ranks -> $SPLIT"
  ibrun -n 8 $FARMS -i elasticity_absorbing.i $COMMON --split-mesh $NR --split-file $SPLIT > $OUT/split.log 2>&1
  echo "[$(date)] split exit $?"
fi
echo "[$(date)] launching $CASE on $SLURM_NNODES nodes, $NR ranks ${REC} EXTRA=[${EXTRA:-}] moose $(git -C $MOOSE_DIR rev-parse --short HEAD)"
ibrun $FARMS -i elasticity_absorbing.i $COMMON --use-split --split-file $SPLIT Outputs/file_base=$OUT/out $REC \
      >> $OUT/run.log 2>&1
rc=$?
echo "[$(date)] $CASE finished (or stopped), farms exit $rc"
exit $rc
