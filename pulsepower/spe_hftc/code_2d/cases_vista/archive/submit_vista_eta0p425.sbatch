#!/bin/bash
# =============================================================================
# SPE-HFTC 2D section on TACC Vista: one 20-pulse train with homogeneous Gc at eta = 0.425 (pulse pressure about
# 210 MPa, between 206 MPa at eta 0.4 and 214 MPa at eta 0.45), absorbing side faces (elasticity_absorbing.i, Lysmer
# dashpots), zeta_pulse = 1e-8 s, quench after 25 us. 72 ranks on one Grace node, the same decomposition as the
# eta 0.4 / 0.45 / 0.5 cases, so that the symmetry breaking of the crack star is seeded the same way.
#   att_homog_eta0p425_absorb   homogeneous Gc (fracture.i, 2721 J/m^2)
# Near-threshold test: eta 0.4 keeps a symmetric six-armed star that arrests at 10.6 mm, eta 0.45 localizes into two
# arms that reach the face (20.0 mm at pulse 20). Does 0.425 stall and then localize (the measured fast-slow-fast
# history) or arrest?
# Submit from this folder:  sbatch submit_vista_eta0p425.sbatch     Continue:  sbatch --export=ALL,RECOVER=1 submit_vista_eta0p425.sbatch
# =============================================================================
#SBATCH -J spe2d_eta0p425
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -t 12:00:00
#SBATCH -A EAR20006
#SBATCH -o spe2d_eta0p425.o%j
#SBATCH -e spe2d_eta0p425.e%j

source $WORK/build_logs/env_moose.sh
FARMS=$FARMS_DIR/farms-opt
CODE=$FARMS_DIR/pulsepower/spe_hftc/code_2d
cd $CODE || exit 1
COMMON="mesh_file=../mesh/spe_section_h0p3.msh eta=0.425 zeta_pulse=1e-8 Np=20"
REC=""; [[ "${RECOVER:-0}" == "1" ]] && REC="--recover"
ONLY="${ONLY:-}"   # kept for the watcher's --recover resubmission (single case anyway)
run() {   # name, input file, fracture sub-app input, rank offset
    local name=$1 input=$2 frac=$3 off=$4
    [[ -n "$ONLY" && "$ONLY" != "$name" ]] && { echo "skipping $name (ONLY=$ONLY)"; return; }
    mkdir -p cases_vista/$name
    echo "[$(date)] launching $name ($input, $frac, eta=0.425) at offset $off"
    ibrun -n 72 -o $off $FARMS -i $input $COMMON MultiApps/fracture/input_files=$frac Outputs/file_base=cases_vista/$name/out $REC > cases_vista/$name/run.log 2>&1 &
}
run att_homog_eta0p425_absorb elasticity_absorbing.i fracture.i 0
wait
echo "[$(date)] case finished (or stopped)"
