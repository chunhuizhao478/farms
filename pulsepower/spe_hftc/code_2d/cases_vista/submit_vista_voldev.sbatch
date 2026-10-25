#!/bin/bash
# NOTE 2026-10-25: elasticity.i and fracture_gaussian.i were moved to ../archive (code_2d keeps elasticity_absorbing.i +
# fracture.i). This script is the record of a finished run; to rerun it, copy the archived input(s) back to code_2d first.
# =============================================================================
# SPE-HFTC 2D section on TACC Vista: the production design with the VOLUMETRIC-DEVIATORIC split
# (Materials/elasticity/decomposition=VOLDEV, Amor: shear energy always drives damage, only volumetric
# compression is protected) instead of the spectral split of att_eta0p5_absorb. Everything else identical:
# Gaussian-scattered Gc (fracture_gaussian.i, cv 12 %), eta = 0.5, absorbing side faces, zeta_pulse = 1e-8 s,
# quench after 25 us, 20 pulses. One case, 72 ranks (half a Grace node, same layout as the production case).
# Submit from this folder:  sbatch submit_vista_voldev.sbatch     Continue:  sbatch --export=ALL,RECOVER=1 submit_vista_voldev.sbatch
# =============================================================================
#SBATCH -J spe2d_voldev
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 72
#SBATCH -t 12:00:00
#SBATCH -A EAR20006
#SBATCH -o spe2d_voldev.o%j
#SBATCH -e spe2d_voldev.e%j

source $WORK/build_logs/env_moose.sh
FARMS=$FARMS_DIR/farms-opt
CODE=$FARMS_DIR/pulsepower/spe_hftc/code_2d
cd $CODE || exit 1
COMMON="mesh_file=../mesh/spe_section_h0p3.msh eta=0.5 zeta_pulse=1e-8 Np=20 MultiApps/fracture/input_files=fracture_gaussian.i Materials/elasticity/decomposition=VOLDEV"
REC=""; [[ "${RECOVER:-0}" == "1" ]] && REC="--recover"
name=att_eta0p5_voldev_absorb
mkdir -p cases_vista/$name
echo "[$(date)] launching $name (elasticity_absorbing.i, VOLDEV split)"
ibrun -n 72 -o 0 $FARMS -i elasticity_absorbing.i $COMMON Outputs/file_base=cases_vista/$name/out $REC > cases_vista/$name/run.log 2>&1
echo "[$(date)] $name finished (or stopped), exit $?"
