#!/bin/bash
# =============================================================================
# SPE-HFTC 2D section on TACC Vista: two 20-pulse trains on one Grace node (144 cores, 72 ranks each),
# both with HOMOGENEOUS toughness (fracture.i, constant Gc = 2721 J/m^2), absorbing side faces
# (elasticity_absorbing.i, Lysmer dashpots), attenuation zeta_pulse = 1e-8 s, quench after 25 us:
#   [1] att_homog_eta0p4_absorb   eta = 0.4
#   [2] att_homog_eta0p5_absorb   eta = 0.5
# Twins of the Gaussian-scatter production case att_eta0p5_absorb (submit_vista_att.sbatch).
# Submit from this folder:  sbatch submit_vista_homog.sbatch     Continue:  sbatch --export=ALL,RECOVER=1 submit_vista_homog.sbatch
# One case only:  sbatch --export=ALL,RECOVER=1,ONLY=att_homog_eta0p4_absorb submit_vista_homog.sbatch
# =============================================================================
#SBATCH -J spe2d_homog
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -t 12:00:00
#SBATCH -A EAR20006
#SBATCH -o spe2d_homog.o%j
#SBATCH -e spe2d_homog.e%j

source $WORK/build_logs/env_moose.sh
FARMS=$FARMS_DIR/farms-opt
CODE=$FARMS_DIR/pulsepower/spe_hftc/code_2d
cd $CODE || exit 1
COMMON="mesh_file=../mesh/spe_section_h0p3.msh zeta_pulse=1e-8 Np=20 MultiApps/fracture/input_files=fracture.i"
REC=""; [[ "${RECOVER:-0}" == "1" ]] && REC="--recover"
ONLY="${ONLY:-}"   # set ONLY=<case name> to launch a single case (e.g. a --recover continuation of one case)
run() {   # name, input file, eta, rank offset
    local name=$1 input=$2 eta=$3 off=$4
    [[ -n "$ONLY" && "$ONLY" != "$name" ]] && { echo "skipping $name (ONLY=$ONLY)"; return; }
    mkdir -p cases_vista/$name
    echo "[$(date)] launching $name ($input, eta=$eta) at offset $off"
    ibrun -n 72 -o $off $FARMS -i $input $COMMON eta=$eta Outputs/file_base=cases_vista/$name/out $REC > cases_vista/$name/run.log 2>&1 &
}
run att_homog_eta0p4_absorb elasticity_absorbing.i 0.4  0
run att_homog_eta0p5_absorb elasticity_absorbing.i 0.5 72
wait
echo "[$(date)] both cases finished (or stopped)"
