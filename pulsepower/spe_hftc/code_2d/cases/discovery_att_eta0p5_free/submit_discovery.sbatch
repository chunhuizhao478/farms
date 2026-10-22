#!/bin/bash
# SPE-HFTC 2D section, Gaussian-scattered toughness, eta = 0.5, attenuation zeta_pulse = 1e-8 s,
# side faces: free (elasticity.i), full 20-pulse train. USC Discovery, oneweek partition (20-core nodes).
# Submit from this folder:  sbatch submit_discovery.sbatch     Continue:  sbatch --export=ALL,RECOVER=1 submit_discovery.sbatch
#SBATCH --job-name=spe2d_att_free
#SBATCH --account=scec_608
#SBATCH --partition=oneweek
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1
#SBATCH --mem=58G
#SBATCH --time=48:00:00
#SBATCH --output=%x.o%j
#SBATCH --error=%x.e%j

module purge
module load usc/13.3.0 gcc/13.3.0 openmpi/5.0.5 openblas/0.3.28 python/3.11.9
FARMS=/home1/zhaochun/projects/farms_cdms/farms-opt
CODE_DIR="$(cd "$SLURM_SUBMIT_DIR/../.." && pwd)"
ARGS="mesh_file=../mesh/spe_section_h0p3.msh eta=0.5 zeta_pulse=1e-8 Np=20 MultiApps/fracture/input_files=fracture_gaussian.i Outputs/file_base=cases/discovery_att_eta0p5_free/out"
cd "$CODE_DIR" || exit 1
echo "job $SLURM_JOB_ID on $(hostname), $SLURM_NTASKS tasks, started $(date)"
if [[ "${RECOVER:-0}" == "1" ]]; then mpirun -np $SLURM_NTASKS "$FARMS" -i elasticity.i $ARGS --recover; else mpirun -np $SLURM_NTASKS "$FARMS" -i elasticity.i $ARGS; fi
echo "finished $(date)"
