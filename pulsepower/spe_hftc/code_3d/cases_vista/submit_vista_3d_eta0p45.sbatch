#!/bin/bash
# =============================================================================
# SPE-HFTC 3D cube on TACC Vista: 20-pulse train, homogeneous Gc, eta = 0.45, absorbing side faces (Lysmer dashpots on
# xmin/xmax/ymin/ymax, holder u_z = 0 on zmin/zmax), production mesh spe_cube_h0p3.msh (14.6 M tets). The 3D twin of the 2D
# section case att_homog_eta0p45_absorb (same inputs apart from the dimension: ../elasticity_absorbing.i + ../fracture.i).
# The run is one continuous train; a job that hits the wall is continued with RECOVER=1 (--recover from the newest
# checkpoint, 100-step interval, 2 kept). Node count and wall time to be confirmed by submit_vista_3d_timing.sbatch first.
# Submit from this folder:  sbatch submit_vista_3d_eta0p45.sbatch      Continue:  sbatch --export=ALL,RECOVER=1 submit_vista_3d_eta0p45.sbatch
# =============================================================================
#SBATCH -J spe3d_eta0p45
#SBATCH -p gg
#SBATCH -N 8
#SBATCH -n 1152
#SBATCH -t 24:00:00
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
COMMON="mesh_file=$MESH eta=0.45 zeta_pulse=1e-8 Np=20"
CASE=cube_homog_eta0p45_absorb
REC=""; [[ "${RECOVER:-0}" == "1" ]] && REC="--recover"
mkdir -p cases_vista/$CASE
if [ ! -d $SPLIT ]; then
  echo "[$(date)] splitting $MESH for $NR ranks -> $SPLIT"
  ibrun -n 8 $FARMS -i elasticity_absorbing.i $COMMON --split-mesh $NR --split-file $SPLIT > cases_vista/$CASE/split.log 2>&1
  echo "[$(date)] split exit $?"
fi
echo "[$(date)] launching $CASE on $SLURM_NNODES nodes, $NR ranks ${REC}"
ibrun $FARMS -i elasticity_absorbing.i $COMMON --use-split --split-file $SPLIT Outputs/file_base=cases_vista/$CASE/out $REC \
      >> cases_vista/$CASE/run.log 2>&1
echo "[$(date)] $CASE finished (or stopped), exit $?"
