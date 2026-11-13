#!/bin/bash
# Stage 1 (LS6 login node, needs internet): clone MOOSE at 6e0969c3, fetch the submodules, build PETSc.
# Run from the farms checkout:  nohup pulsepower/spe_hftc/code_3d/cases_ls6/build_stage1_login.sh > $WORK/spe3d/build_logs/stage1.log 2>&1 &
set -o pipefail
HERE=$(cd $(dirname $0) && pwd)
source $HERE/env_ls6.sh
L=$SPE3D_ROOT/build_logs; mkdir -p $L
echo "[$(date)] modules: $(module -t list 2>&1 | tr '\n' ' ')"; which mpicc; mpicc -show
MOOSE_COMMIT=6e0969c339df46d6815fbb7cf36f0c26e8a46808
if [ ! -d $MOOSE_DIR/.git ]; then
  git clone -q https://github.com/idaholab/moose.git $MOOSE_DIR || exit 2
fi
cd $MOOSE_DIR && git checkout -q $MOOSE_COMMIT || exit 2
echo "[$(date)] moose $(git rev-parse HEAD)"
git submodule update --init --recursive petsc libmesh framework/contrib/wasp > $L/submodules.log 2>&1 || { echo "ERROR submodules"; tail $L/submodules.log; exit 2; }
for s in petsc libmesh framework/contrib/wasp; do echo "$s $(git -C $s rev-parse HEAD)"; done
echo "[$(date)] PETSc build start"
export MOOSE_JOBS=4
./scripts/update_and_rebuild_petsc.sh --skip-submodule-update > $L/petsc.log 2>&1
echo "[$(date)] PETSc exit $?"; tail -5 $L/petsc.log
ls -la petsc/arch-moose/lib/libpetsc.so* || exit 3
ldd petsc/arch-moose/lib/libpetsc.so | grep -E "not found|libmpi"
echo "[$(date)] stage 1 done; next: sbatch build_stage2_libmesh.sbatch"
