#!/bin/bash
# =============================================================================
# SPE-HFTC 2D section on TACC Vista: two 20-pulse trains on one Grace node (144 cores, 72 ranks each):
#   [1] free side faces          elasticity.i
#   [2] absorbing side faces     elasticity_absorbing.i  (Lysmer dashpots on xmin/xmax/ymin/ymax)
# both with Gaussian-scattered toughness, eta = 0.5, attenuation zeta_pulse = 1e-8 s, quench after 25 us.
# Submit from this folder:  sbatch submit_vista_att.sbatch     Continue:  sbatch --export=ALL,RECOVER=1 submit_vista_att.sbatch
# One case only (e.g. continue the absorbing case alone):  sbatch --export=ALL,RECOVER=1,ONLY=att_eta0p5_absorb submit_vista_att.sbatch
# =============================================================================
#SBATCH -J spe2d_att
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -t 12:00:00
#SBATCH -A EAR20006
#SBATCH -o spe2d_att.o%j
#SBATCH -e spe2d_att.e%j

source $WORK/build_logs/env_moose.sh
FARMS=$FARMS_DIR/farms-opt
CODE=$FARMS_DIR/pulsepower/spe_hftc/code_2d
cd $CODE || exit 1
COMMON="mesh_file=../mesh/spe_section_h0p3.msh eta=0.5 zeta_pulse=1e-8 Np=20 MultiApps/fracture/input_files=fracture_gaussian.i"
REC=""; [[ "${RECOVER:-0}" == "1" ]] && REC="--recover"
ONLY="${ONLY:-}"   # set ONLY=<case name> to launch a single case (e.g. a --recover continuation of one case)
run() {   # name, input file, rank offset
    local name=$1 input=$2 off=$3
    [[ -n "$ONLY" && "$ONLY" != "$name" ]] && { echo "skipping $name (ONLY=$ONLY)"; return; }
    mkdir -p cases_vista/$name
    echo "[$(date)] launching $name ($input) at offset $off"
    ibrun -n 72 -o $off $FARMS -i $input $COMMON Outputs/file_base=cases_vista/$name/out $REC > cases_vista/$name/run.log 2>&1 &
}
run att_eta0p5_free   elasticity.i            0
run att_eta0p5_absorb elasticity_absorbing.i 72
wait
echo "[$(date)] both cases finished (or stopped)"
