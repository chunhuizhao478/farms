#!/bin/bash
# =============================================================================
# SPE-HFTC 3D cube on TACC Lonestar6: the Vista case cube_homog_eta0p45_absorb run again from t = 0 (20-pulse train,
# homogeneous Gc, eta = 0.45, absorbing side faces, mesh spe_cube_h0p3.msh, 14.6 M tets; see ../cases_vista/README.md).
# Same inputs and overrides as ../cases_vista/submit_vista_3d_eta0p45.sbatch (dt_quiet = 5e-7 through EXTRA).
# 18 nodes x 128 cores = 2304 ranks, the rank count of the Vista run, so the existing split spe_cube_h0p3_split2304.cpr
# (linked by prepare_ls6_case.sh) is reused. Outputs, run.log and checkpoints (about 54 GB each) go to
# $SCRATCH/spe3d/<case>, never to the checkout: /work has a 1 TB quota shared with Vista and filled up on 2026-11-15.
# Checkpoints every 100 steps, 2 kept; a job that hits the 48 h wall is
# continued with RECOVER=1. The job exits with the farms exit code, so an afternotok continuation also starts after a crash.
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
MESH=../mesh/spe_cube_h0p3.msh
SPLIT=../mesh/spe_cube_h0p3_split${NR}.cpr
COMMON="mesh_file=$MESH eta=0.45 zeta_pulse=1e-8 Np=20 ${EXTRA:-}"
CASE=cube_homog_eta0p45_absorb
REC=""; [[ "${RECOVER:-0}" == "1" ]] && REC="--recover"
OUT=$SCRATCH/spe3d/$CASE
mkdir -p $OUT || exit 2
[ -e $SPLIT ] || { echo "ERROR: $SPLIT missing (run prepare_ls6_case.sh)"; exit 2; }
echo "[$(date)] launching $CASE on $SLURM_NNODES nodes, $NR ranks ${REC} EXTRA=[${EXTRA:-}] moose $(git -C $MOOSE_DIR rev-parse --short HEAD)"
ibrun $FARMS -i elasticity_absorbing.i $COMMON --use-split --split-file $SPLIT Outputs/file_base=$OUT/out $REC \
      >> $OUT/run.log 2>&1
rc=$?
echo "[$(date)] $CASE finished (or stopped), farms exit $rc"
exit $rc
