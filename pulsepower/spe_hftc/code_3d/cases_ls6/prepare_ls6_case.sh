#!/bin/bash
# Stage the files of the 3D cube run that are not in git, from the Vista tree on the shared /work file system:
# the mesh, the 2304-rank split, and checkpoint 1370 of the Vista run (job 1104488, t = 239.4 us, pulse 4 of 20).
# Mesh and split are symlinks. The checkpoint is hard-linked (no extra space, same file system), because the Vista run
# keeps only 2 checkpoints and deletes 1370 once it has advanced; a hard link survives that. Run on an LS6 login node.
HERE=$(cd $(dirname $0) && pwd)
source $HERE/env_ls6.sh
V=/work/10024/zhaochun/vista/farms_cdms/pulsepower/spe_hftc
C=$FARMS_DIR/pulsepower/spe_hftc
CASE=$C/code_3d/cases_ls6/cube_homog_eta0p45_absorb
STEP=${STEP:-1370}
mkdir -p $C/mesh $CASE/out_cp
for f in spe_cube_h0p3.msh spe_cube_h0p3_split2304.cpr; do
  [ -e $V/mesh/$f ] || { echo "ERROR: $V/mesh/$f missing"; exit 2; }
  ln -sfn $V/mesh/$f $C/mesh/$f
done
SRC=$V/code_3d/cases_vista/cube_homog_eta0p45_absorb/out_cp
ls -d $SRC/$STEP-mesh.cpa.gz > /dev/null || { echo "ERROR: checkpoint $STEP missing in $SRC"; exit 2; }
if [ ! -e $CASE/out_cp/$STEP-mesh.cpa.gz ]; then
  cp -al $SRC/$STEP-mesh.cpa.gz $CASE/out_cp/ || exit 3
  for d in $SRC/$STEP-restart-*.rd; do cp -al $d $CASE/out_cp/ || exit 3; done
fi
echo "mesh: $(ls -la $C/mesh | grep -c h0p3) links; checkpoint $STEP: $(ls -d $CASE/out_cp/$STEP-restart-*.rd | wc -l) rank folders (expect 2304)"
