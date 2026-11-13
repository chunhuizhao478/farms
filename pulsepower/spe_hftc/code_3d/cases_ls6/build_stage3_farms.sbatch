#!/bin/bash
# Stage 3 (LS6 compute node): MOOSE framework + modules + farms.
#SBATCH -J spe3d_farms
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH -A EAR24034
#SBATCH -o farms.o%j
#SBATCH -e farms.e%j
source $WORK/spe3d/farms_cdms/pulsepower/spe_hftc/code_3d/cases_ls6/env_ls6.sh
source $FARMS_DIR/pulsepower/spe_hftc/code_3d/cases_ls6/guard.sh
export MOOSE_JOBS=64
cd $FARMS_DIR
echo "=== farms (framework + modules + app) with MOOSE_DIR=$MOOSE_DIR ($(git -C $MOOSE_DIR rev-parse --short HEAD))  $(date) ==="
make -j 64 > $SPE3D_ROOT/build_logs/farms_make.log 2>&1
tail -15 $SPE3D_ROOT/build_logs/farms_make.log
ls -la farms-opt && ./farms-opt --version | tail -1 || { echo "FARMS BUILD FAILED"; exit 1; }
ldd farms-opt | grep -E "not found|libmpi"
echo "=== done $(date), farms commit $(git rev-parse HEAD) ==="
