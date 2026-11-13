#!/bin/bash
# Stage the mesh of the 3D cube run (not in git) from the Vista tree on the shared /work file system: the .msh file and
# the 2304-rank split, as symlinks. The LS6 run starts from t = 0 (decision 2026-11-13), so no Vista checkpoint is staged.
# Run on an LS6 login node.
HERE=$(cd $(dirname $0) && pwd)
source $HERE/env_ls6.sh
V=/work/10024/zhaochun/vista/farms_cdms/pulsepower/spe_hftc
C=$FARMS_DIR/pulsepower/spe_hftc
mkdir -p $C/mesh
for f in spe_cube_h0p3.msh spe_cube_h0p3_split2304.cpr; do
  [ -e $V/mesh/$f ] || { echo "ERROR: $V/mesh/$f missing"; exit 2; }
  ln -sfn $V/mesh/$f $C/mesh/$f
done
ls -la $C/mesh | grep h0p3
