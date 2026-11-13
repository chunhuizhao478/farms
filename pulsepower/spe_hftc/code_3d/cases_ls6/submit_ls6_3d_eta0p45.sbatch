#!/bin/bash
# =============================================================================
# SPE-HFTC 3D cube on TACC Lonestar6: continuation of the Vista run cube_homog_eta0p45_absorb (20-pulse train,
# homogeneous Gc, eta = 0.45, absorbing side faces, mesh spe_cube_h0p3.msh, 14.6 M tets; see ../cases_vista/README.md).
# Vista segment 1 (job 1104488) reached step 1371, t = 239.9 us (pulse 4 of 20); its checkpoint 1370 is staged by
# prepare_ls6_case.sh. A recover needs the rank count of the checkpoint: 2304 ranks = 18 LS6 nodes x 128 cores.
# Inputs and overrides are those of ../cases_vista/submit_vista_3d_eta0p45.sbatch (dt_quiet = 5e-7 through EXTRA).
# The job exits with the farms exit code, so that an afternotok continuation also starts after a crash.
# Submit from this folder (after prepare_ls6_case.sh; RECOVER=1 picks the newest checkpoint in the case folder):
#   sbatch --export=ALL,RECOVER=1,EXTRA="dt_quiet=5e-7" submit_ls6_3d_eta0p45.sbatch
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
mkdir -p cases_ls6/$CASE
[ -e $SPLIT ] || { echo "ERROR: $SPLIT missing (run prepare_ls6_case.sh)"; exit 2; }
echo "[$(date)] launching $CASE on $SLURM_NNODES nodes, $NR ranks ${REC} EXTRA=[${EXTRA:-}] moose $(git -C $MOOSE_DIR rev-parse --short HEAD)"
ibrun $FARMS -i elasticity_absorbing.i $COMMON --use-split --split-file $SPLIT Outputs/file_base=cases_ls6/$CASE/out $REC \
      >> cases_ls6/$CASE/run.log 2>&1
rc=$?
echo "[$(date)] $CASE finished (or stopped), farms exit $rc"
exit $rc
