#!/bin/bash
# =============================================================================
# SPE-HFTC 3D cube on TACC Vista: 20-pulse train, homogeneous Gc, eta = 0.45, absorbing side faces (Lysmer dashpots on
# xmin/xmax/ymin/ymax, holder u_z = 0 on zmin/zmax), production mesh spe_cube_h0p3.msh (14.6 M tets). The 3D twin of the 2D
# section case att_homog_eta0p45_absorb (same inputs apart from the dimension: ../elasticity_absorbing.i + ../fracture.i).
# The run is one continuous train; a job that hits the wall is continued with RECOVER=1 (--recover from the newest
# checkpoint, 100-step interval, 2 kept). Node count from the timing runs 1093797/1093056/1093798 (65/42/33 s per step on
# 4/8/16 nodes, 77 % then 64 % efficiency): 16 nodes. EXTRA is appended to the farms command line (input overrides), e.g.
# EXTRA="dt_quiet=5e-7" for the coarser ring-down step (see cases_vista/README.md); the input defaults are unchanged.
# Submit from this folder:  sbatch --export=ALL,EXTRA="dt_quiet=5e-7" submit_vista_3d_eta0p45.sbatch
# Continue:                 sbatch --export=ALL,RECOVER=1,EXTRA="dt_quiet=5e-7" submit_vista_3d_eta0p45.sbatch
# Another node count: add -N NN -n $((NN*144)) (the split ../mesh/spe_cube_h0p3_split<ranks>.cpr is built if missing).
# =============================================================================
#SBATCH -J spe3d_eta0p45
#SBATCH -p gg
#SBATCH -N 16
#SBATCH -n 2304
#SBATCH -t 48:00:00
#SBATCH -A EAR20006
#SBATCH -o spe3d_eta0p45.o%j
#SBATCH -e spe3d_eta0p45.e%j

source $WORK/build_logs/env_moose.sh
FARMS=$FARMS_DIR/farms-opt
CODE=$FARMS_DIR/pulsepower/spe_hftc/code_3d
cd $CODE || exit 1
NR=$SLURM_NTASKS
MESH=../mesh/spe_cube_h0p3.msh
SPLIT=../mesh/spe_cube_h0p3_split${NR}.cpr
COMMON="mesh_file=$MESH eta=0.45 zeta_pulse=1e-8 Np=20 ${EXTRA:-}"
CASE=cube_homog_eta0p45_absorb
REC=""; [[ "${RECOVER:-0}" == "1" ]] && REC="--recover"
mkdir -p cases_vista/$CASE
if [ ! -d $SPLIT ]; then
  echo "[$(date)] splitting $MESH for $NR ranks -> $SPLIT"
  ibrun -n 8 $FARMS -i elasticity_absorbing.i $COMMON --split-mesh $NR --split-file $SPLIT > cases_vista/$CASE/split.log 2>&1
  echo "[$(date)] split exit $?"
fi
echo "[$(date)] launching $CASE on $SLURM_NNODES nodes, $NR ranks ${REC} EXTRA=[${EXTRA:-}]"
ibrun $FARMS -i elasticity_absorbing.i $COMMON --use-split --split-file $SPLIT Outputs/file_base=cases_vista/$CASE/out $REC \
      >> cases_vista/$CASE/run.log 2>&1
echo "[$(date)] $CASE finished (or stopped), exit $?"
