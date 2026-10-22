#!/bin/bash
# =============================================================================
# SPE-HFTC 2D section on TACC Vista (Grace CPU node, 144 cores): four cases at once,
# 36 ranks each, full 20-pulse train. Gaussian-scattered toughness; attenuation
# zeta_pulse = 1e-8 s on three of them (the face-spall remedy), one without for reference.
# Submit from this folder:          sbatch submit_vista_sweep.sbatch
# Continue after a time-out:        sbatch --export=ALL,RECOVER=1 submit_vista_sweep.sbatch
# =============================================================================
#SBATCH -J spe2d_sweep
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -t 12:00:00
#SBATCH -A EAR20006
#SBATCH -o spe2d_sweep.o%j
#SBATCH -e spe2d_sweep.e%j

source $WORK/build_logs/env_moose.sh                 # gcc 14.2 + gcc14 OpenMPI 5.0.5, MOOSE_DIR, FARMS_DIR
FARMS=$FARMS_DIR/farms-opt
CODE=$FARMS_DIR/pulsepower/spe_hftc/code_2d
CASES=$CODE/cases_vista
cd $CODE || exit 1
COMMON="mesh_file=../mesh/spe_section_h0p3.msh Np=20 MultiApps/fracture/input_files=fracture_gaussian.i"
REC=""; [[ "${RECOVER:-0}" == "1" ]] && REC="--recover"

run() {   # name, rank offset, overrides...
    local name=$1 off=$2; shift 2
    mkdir -p $CASES/$name
    echo "[$(date)] launching $name at offset $off: $*"
    ibrun -n 36 -o $off $FARMS -i elasticity.i $COMMON "$@" Outputs/file_base=cases_vista/$name/out $REC > $CASES/$name/run.log 2>&1 &
}
run eta0p4_att     0 eta=0.4 zeta_pulse=1e-8
run eta0p5_att    36 eta=0.5 zeta_pulse=1e-8
run eta0p6_att    72 eta=0.6 zeta_pulse=1e-8
run eta0p5_noatt 108 eta=0.5
wait
echo "[$(date)] all four cases finished (or stopped)"
