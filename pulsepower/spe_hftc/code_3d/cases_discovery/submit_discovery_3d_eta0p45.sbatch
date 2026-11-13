#!/bin/bash
# =============================================================================
# SPE-HFTC 3D cube on USC Discovery: the Vista case cube_homog_eta0p45_absorb run from t = 0 (20-pulse train,
# homogeneous Gc, eta = 0.45, absorbing side faces, mesh spe_cube_h0p3.msh, 14.6 M tets; see ../cases_vista/README.md).
# Same inputs and overrides as ../cases_vista/submit_vista_3d_eta0p45.sbatch (dt_quiet = 5e-7 through EXTRA).
# 16 EPYC nodes x 64 cores = 1024 ranks (the normal QOS caps a user at 2000 CPUs, so the 2304-rank Vista checkpoints
# cannot be recovered here). The mesh split for the rank count is built by the job if missing. Outputs, checkpoints
# (100-step interval, 2 kept) and the split go to /scratch1 (they do not fit the /home1 quota); the mesh file
# ../mesh/spe_cube_h0p3.msh is not in git and must be copied into the checkout first.
# The job exits with the farms exit code, so an afternotok continuation also starts after a crash.
# Submit from this folder:  sbatch --export=ALL,EXTRA="dt_quiet=5e-7" submit_discovery_3d_eta0p45.sbatch
# Continue:                 sbatch --export=ALL,RECOVER=1,EXTRA="dt_quiet=5e-7" submit_discovery_3d_eta0p45.sbatch
# =============================================================================
#SBATCH --job-name=spe3d_eta0p45
#SBATCH --account=scec_608
#SBATCH --partition=main
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --time=48:00:00
#SBATCH --output=%x.o%j
#SBATCH --error=%x.e%j

module purge
module load usc/13.3.0 gcc/13.3.0 openmpi/5.0.5 openblas/0.3.28 python/3.11.9
FARMS_DIR=/home1/zhaochun/projects/farms_cdms
cd $FARMS_DIR || exit 2
git fetch -q origin || { echo "ERROR: git fetch failed"; exit 2; }
[ -z "$(git status --porcelain --untracked-files=no)" ] || { echo "ERROR: uncommitted changes in $FARMS_DIR"; git status --short; exit 2; }
[ -n "$(git branch -r --contains HEAD)" ] || { echo "ERROR: HEAD $(git rev-parse --short HEAD) is not on GitHub"; exit 2; }
echo "farms commit $(git rev-parse HEAD)"
FARMS=$FARMS_DIR/farms-opt
cd $FARMS_DIR/pulsepower/spe_hftc/code_3d || exit 1
NR=$SLURM_NTASKS
CASE=cube_homog_eta0p45_absorb
OUT=/scratch1/$USER/spe_hftc_3d
MESH=../mesh/spe_cube_h0p3.msh
SPLIT=$OUT/spe_cube_h0p3_split${NR}.cpr
COMMON="mesh_file=$MESH eta=0.45 zeta_pulse=1e-8 Np=20 ${EXTRA:-}"
REC=""; [[ "${RECOVER:-0}" == "1" ]] && REC="--recover"
mkdir -p $OUT/$CASE
[ -e $MESH ] || { echo "ERROR: $MESH missing"; exit 2; }
if [ ! -d $SPLIT ]; then
  echo "[$(date)] splitting $MESH for $NR ranks -> $SPLIT"
  mpirun -np 8 $FARMS -i elasticity_absorbing.i $COMMON --split-mesh $NR --split-file $SPLIT > $OUT/$CASE/split.log 2>&1
  echo "[$(date)] split exit $?"
fi
echo "[$(date)] launching $CASE on $SLURM_NNODES nodes, $NR ranks ${REC} EXTRA=[${EXTRA:-}]"
mpirun -np $NR $FARMS -i elasticity_absorbing.i $COMMON --use-split --split-file $SPLIT Outputs/file_base=$OUT/$CASE/out $REC \
      >> $OUT/$CASE/run.log 2>&1
rc=$?
echo "[$(date)] $CASE finished (or stopped), farms exit $rc"
exit $rc
