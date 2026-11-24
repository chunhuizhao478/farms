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
# 2026-11-16: outputs, run.log and checkpoints (about 54 GB each) go to $SCRATCH/spe3d_vista/<case>, never to the checkout
# (/work has a 1 TB quota that filled up on 2026-11-15). The inputs are read from the checkout this script is submitted
# from (clean and on GitHub: guard below); the binary is $FARMS_DIR/farms-opt of env_moose.sh. RECOVER_FROM=<dir>/<step>
# recovers from an explicit checkpoint (e.g. the step-1370 checkpoint of segment 1); RECOVER=1 from the newest one in the
# case folder. The job exits with the farms exit code. run_chain_vista.sh drives the segments with a stall guard.
# 2026-11-17 geometry: 7.94 mm borehole, electrode gap at mid-height, mesh spe_cube_d7p94_h0p3.msh, case
# cube_d7p94_mid_eta0p45_absorb from t = 0 (the header above describes the 8.3 mm / Lz/3 runs up to 2026-11-16).
# The split is built by the job on scratch if missing.
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
CODE=$(cd $SLURM_SUBMIT_DIR/.. && pwd)
cd $CODE || exit 1
git fetch -q origin || { echo "ERROR: git fetch failed"; exit 2; }
[ -z "$(git status --porcelain --untracked-files=no)" ] || { echo "ERROR: uncommitted changes in $CODE"; git status --short; exit 2; }
[ -n "$(git branch -r --contains HEAD)" ] || { echo "ERROR: HEAD $(git rev-parse --short HEAD) is not on GitHub"; exit 2; }
echo "farms inputs commit $(git rev-parse HEAD), binary $FARMS"
NR=$SLURM_NTASKS
MESH=../mesh/spe_cube_d7p94_h0p3.msh
SPLIT=$SCRATCH/spe3d_vista/spe_cube_d7p94_h0p3_split${NR}.cpr
COMMON="mesh_file=$MESH eta=0.45 zeta_pulse=1e-8 Np=20 ${EXTRA:-}"
CASE=cube_d7p94_mid_eta0p45_absorb
OUT=$SCRATCH/spe3d_vista/$CASE
REC=""; [[ "${RECOVER:-0}" == "1" ]] && REC="--recover"
[ -n "${RECOVER_FROM:-}" ] && REC="--recover $RECOVER_FROM"
mkdir -p $OUT || exit 2
[ -e $MESH ] || { echo "ERROR: $MESH missing (copy it into the checkout, it is not in git)"; exit 2; }
if [ ! -d $SPLIT ]; then
  echo "[$(date)] splitting $MESH for $NR ranks -> $SPLIT"
  ibrun -n 8 $FARMS -i elasticity_absorbing.i $COMMON --split-mesh $NR --split-file $SPLIT > $OUT/split.log 2>&1
  echo "[$(date)] split exit $?"
fi
# 2026-11-24 node-local staging: the stalls of jobs 1159778 (i615-004) and 1164089 (i617-063) were one rank in D state
# (folio_wait_bit_common / ptlrpc_set_wait = a Lustre page-in) while all others waited in VecAssemblyBegin: code pages of
# farms-opt's libraries on /work (Lustre) were re-read mid-run and the read hung. sbcast copies farms-opt and its shared
# libraries (Slurm skips /lib,/usr/lib,/lib64,/usr/lib64) to /tmp on every node; a wrapper puts them first on
# LD_LIBRARY_PATH (farms-opt and the MOOSE libraries use RUNPATH, which LD_LIBRARY_PATH overrides). Falls back to the
# /work binary if staging fails. STAGE=0 turns it off.
LAUNCH=$FARMS
if [ "${STAGE:-1}" = 1 ]; then
  SEXE=/tmp/farms-opt;   # keep the name farms-opt: the stall watchdog finds ranks with pgrep -x farms-opt
  SWRAP=/tmp/farms-run.$SLURM_JOB_ID.sh
  printf '#!/bin/bash\nexport LD_LIBRARY_PATH=%s_libs:$LD_LIBRARY_PATH\nexec %s "$@"\n' $SEXE $SEXE > $SWRAP; chmod 755 $SWRAP
  if sbcast --force --send-libs=yes $FARMS $SEXE && sbcast --force $SWRAP $SWRAP; then
    nw=$(LD_LIBRARY_PATH=${SEXE}_libs:$LD_LIBRARY_PATH ldd $SEXE | grep -c '=> /work')
    echo "[$(date)] staged farms-opt + $(ls ${SEXE}_libs | wc -l) libraries to /tmp on $SLURM_NNODES nodes ($(du -sh ${SEXE}_libs | cut -f1)); libraries still resolved from /work: $nw"
    LAUNCH=$SWRAP
  else
    echo "[$(date)] WARNING: sbcast staging failed, running farms-opt from /work"
  fi
fi
echo "[$(date)] launching $CASE on $SLURM_NNODES nodes, $NR ranks ${REC} EXTRA=[${EXTRA:-}] exe=$LAUNCH"
# 10 min in: count /work-mapped shared objects in one running rank on this node (0 = staging works)
( sleep 600; p=$(pgrep -u $USER -x farms-opt | head -1); [ -n "$p" ] && echo "[$(date)] rank pid $p maps $(grep -c '/work/.*\.so' /proc/$p/maps) /work .so mappings, $(grep -c '/tmp/farms-opt' /proc/$p/maps) /tmp ones" ) &
ibrun $LAUNCH -i elasticity_absorbing.i $COMMON --use-split --split-file $SPLIT Outputs/file_base=$OUT/out $REC \
      >> $OUT/run.log 2>&1
rc=$?
echo "[$(date)] $CASE finished (or stopped), farms exit $rc"
exit $rc
