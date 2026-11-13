#!/bin/bash
# Stage 2 (LS6 compute node, no internet): libMesh + WASP against the PETSc of stage 1.
#SBATCH -J spe3d_libmesh
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH -A EAR24034
#SBATCH -o libmesh.o%j
#SBATCH -e libmesh.e%j
source $WORK/spe3d/farms_cdms/pulsepower/spe_hftc/code_3d/cases_ls6/env_ls6.sh
source $FARMS_DIR/pulsepower/spe_hftc/code_3d/cases_ls6/guard.sh
module -t list 2>&1 | tr '\n' ' '; echo; which mpicxx; mpicxx --version | head -1
export MOOSE_JOBS=64
cd $MOOSE_DIR; echo "moose $(git rev-parse HEAD)"
echo "=== libMesh (methods: $METHODS) $(date) ==="
./scripts/update_and_rebuild_libmesh.sh --skip-submodule-update --disable-netgen > $SPE3D_ROOT/build_logs/libmesh_build_full.log 2>&1
tail -5 $SPE3D_ROOT/build_logs/libmesh_build_full.log
ls -la libmesh/installed/lib/libmesh_opt.so* || { echo "LIBMESH FAILED"; exit 1; }
echo "=== WASP $(date) ==="
./scripts/update_and_rebuild_wasp.sh --skip-submodule-update > $SPE3D_ROOT/build_logs/wasp.log 2>&1
ls -la framework/contrib/wasp/install/lib/libwaspcore.so* || { echo "WASP FAILED"; tail -20 $SPE3D_ROOT/build_logs/wasp.log; exit 1; }
ldd libmesh/installed/lib/libmesh_opt.so | grep -E "not found|libmpi"
echo "=== done $(date) ==="
