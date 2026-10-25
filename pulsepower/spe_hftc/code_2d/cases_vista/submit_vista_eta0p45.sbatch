#!/bin/bash
# NOTE 2026-10-25: elasticity.i and fracture_gaussian.i were moved to ../archive (code_2d keeps elasticity_absorbing.i +
# fracture.i). This script is the record of a finished run; to rerun it, copy the archived input(s) back to code_2d first.
# =============================================================================
# SPE-HFTC 2D section on TACC Vista: two 20-pulse trains on one Grace node (144 cores, 72 ranks each),
# both with eta = 0.45 (pulse pressure 214 MPa against 222.5 MPa at eta 0.5 and 206 MPa at eta 0.4),
# absorbing side faces (elasticity_absorbing.i, Lysmer dashpots), zeta_pulse = 1e-8 s, quench after 25 us:
#   [1] att_eta0p45_absorb         Gaussian Gc scatter (fracture_gaussian.i, cv 0.12, seed 1): the production design at eta 0.45
#   [2] att_homog_eta0p45_absorb   homogeneous Gc (fracture.i, 2721 J/m^2): the twin
# Completes the sweep eta = 0.4 / 0.45 / 0.5 (att_homog_eta0p4_absorb, att_homog_eta0p5_absorb, att_eta0p5_absorb).
# Submit from this folder:  sbatch submit_vista_eta0p45.sbatch     Continue:  sbatch --export=ALL,RECOVER=1 submit_vista_eta0p45.sbatch
# One case only:  sbatch --export=ALL,RECOVER=1,ONLY=att_eta0p45_absorb submit_vista_eta0p45.sbatch
# =============================================================================
#SBATCH -J spe2d_eta0p45
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -t 12:00:00
#SBATCH -A EAR20006
#SBATCH -o spe2d_eta0p45.o%j
#SBATCH -e spe2d_eta0p45.e%j

source $WORK/build_logs/env_moose.sh
FARMS=$FARMS_DIR/farms-opt
CODE=$FARMS_DIR/pulsepower/spe_hftc/code_2d
cd $CODE || exit 1
COMMON="mesh_file=../mesh/spe_section_h0p3.msh eta=0.45 zeta_pulse=1e-8 Np=20"
REC=""; [[ "${RECOVER:-0}" == "1" ]] && REC="--recover"
ONLY="${ONLY:-}"   # set ONLY=<case name> to launch a single case (e.g. a --recover continuation of one case)
run() {   # name, input file, fracture sub-app input, rank offset
    local name=$1 input=$2 frac=$3 off=$4
    [[ -n "$ONLY" && "$ONLY" != "$name" ]] && { echo "skipping $name (ONLY=$ONLY)"; return; }
    mkdir -p cases_vista/$name
    echo "[$(date)] launching $name ($input, $frac, eta=0.45) at offset $off"
    ibrun -n 72 -o $off $FARMS -i $input $COMMON MultiApps/fracture/input_files=$frac Outputs/file_base=cases_vista/$name/out $REC > cases_vista/$name/run.log 2>&1 &
}
run att_eta0p45_absorb       elasticity_absorbing.i fracture_gaussian.i  0
run att_homog_eta0p45_absorb elasticity_absorbing.i fracture.i          72
wait
echo "[$(date)] both cases finished (or stopped)"
