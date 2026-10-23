#!/bin/bash
# SPE-HFTC 2D section on TACC Vista: time-step sensitivity check for the 3D schedule, two 20-pulse trains on one Grace node
# (72 ranks each), both the homogeneous-Gc eta = 0.45 absorbing case (the 3D reference, baseline att_homog_eta0p45_absorb
# with dt_max = dt_quiet = 1e-7, 12.5 k steps):
#   [A] att_homog_eta0p45_absorb_dtq5      dt_quiet = 5e-7 in the damped ring-down (t > 25 us of each window); growth phase unchanged
#   [B] att_homog_eta0p45_absorb_dtq5_dtm2 dt_quiet = 5e-7 and dt_max = 2e-7 in the growth phase (crack tip <= 0.8 h per step)
# Pass criterion: extents at pulses 5/11/15/20 and the crack energy within the one-realization scatter of the baseline.
# Submit from this folder:  sbatch submit_vista_dtcheck.sbatch     Continue:  sbatch --export=ALL,RECOVER=1 submit_vista_dtcheck.sbatch
# One case only:  sbatch --export=ALL,RECOVER=1,ONLY=att_homog_eta0p45_absorb_dtq5 submit_vista_dtcheck.sbatch
#SBATCH -J spe2d_dtcheck
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -t 10:00:00
#SBATCH -A EAR20006
#SBATCH -o spe2d_dtcheck.o%j
#SBATCH -e spe2d_dtcheck.e%j

source $WORK/build_logs/env_moose.sh
FARMS=$FARMS_DIR/farms-opt
CODE=$FARMS_DIR/pulsepower/spe_hftc/code_2d
cd $CODE || exit 1
COMMON="mesh_file=../mesh/spe_section_h0p3.msh eta=0.45 zeta_pulse=1e-8 Np=20"
REC=""; [[ "${RECOVER:-0}" == "1" ]] && REC="--recover"
ONLY="${ONLY:-}"
run() {   # name, extra parameters, rank offset
    local name=$1 extra=$2 off=$3
    [[ -n "$ONLY" && "$ONLY" != "$name" ]] && { echo "skipping $name (ONLY=$ONLY)"; return; }
    mkdir -p cases_vista/$name
    echo "[$(date)] launching $name ($extra) at offset $off"
    ibrun -n 72 -o $off $FARMS -i elasticity_absorbing.i $COMMON $extra MultiApps/fracture/input_files=fracture.i \
          Outputs/file_base=cases_vista/$name/out $REC > cases_vista/$name/run.log 2>&1 &
}
run att_homog_eta0p45_absorb_dtq5      "dt_quiet=5e-7"            0
run att_homog_eta0p45_absorb_dtq5_dtm2 "dt_quiet=5e-7 dt_max=2e-7" 72
wait
echo "[$(date)] both cases finished (or stopped)"
