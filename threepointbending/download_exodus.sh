#!/bin/bash
# Download all exodus files from Frontera

REMOTE="zhaochun@frontera.tacc.utexas.edu"
REMOTE_BASE="/scratch2/10024/zhaochun/projects/farms_cdms_01292026/threepointbending/code_angles"
LOCAL_DIR="/Volumes/One Touch/Research/ThreePointBending"

for label in 0.00 0.25 0.50 0.75 1.00 1.25 1.50 1.75 2.00 2.25 2.50 2.75 3.00; do
  name="case_whole_3d_x${label}mm_y1.6mm"
  echo "Downloading ${name}..."
  scp "${REMOTE}:${REMOTE_BASE}/${name}/elasticity_exodus.e" "${LOCAL_DIR}/${name}.e"
done
