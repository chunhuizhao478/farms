#!/bin/bash
# The mesh of the 3D cube run is not in git: copy ../../mesh/spe_cube_d7p94_h0p3.msh (built with make_spe_cube_mesh.py
# --rb 3.97e-3 ..., see ../../mesh/README.md) into the checkout on LS6, then run this check. The split is built by the job.
HERE=$(cd $(dirname $0) && pwd)
source $HERE/env_ls6.sh
M=$FARMS_DIR/pulsepower/spe_hftc/mesh/spe_cube_d7p94_h0p3.msh
[ -s $M ] || { echo "ERROR: $M missing"; exit 2; }
ls -la $M; head -2 $M | tail -1
