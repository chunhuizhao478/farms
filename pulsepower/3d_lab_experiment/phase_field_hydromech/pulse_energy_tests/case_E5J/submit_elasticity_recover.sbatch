#!/bin/bash
#SBATCH -J lab_3dhm_E5J                        # Job name
#SBATCH -o lab_3dhm_E5J.o%j                    # Name of stdout output file
#SBATCH -e lab_3dhm_E5J.e%j                    # Name of stderr error file
#SBATCH -p normal                              # Queue (partition) name
#SBATCH -N 10                                  # Total # of nodes
#SBATCH -n 500                                 # Total # of mpi tasks
#SBATCH -t 48:00:00                            # Run time (hh:mm:ss)
#SBATCH --mail-type=all                        # Send email at begin and end of job
#SBATCH -A EAR20006                            # Project/Allocation name
#SBATCH --mail-user=chunhui3@illinois.edu


# Load necessary modules
module swap intel gcc
module load cuda
export CXXFLAGS=-I/opt/apps/gcc/9.1.0/include/c++/9.1.0/
export CC=mpicc CXX=mpicxx FC=mpif90 F90=mpif90 F77=mpif77

# Enable MPI debugging
export MV2_DEBUG=1
export MV2_SHOW_ENV_INFO=1

export MOOSE_JOBS=6 METHODS=opt

# Run the simulation
# ibrun -n 1 /scratch2/10024/zhaochun/projects/farms_cdms_04192026/farms-opt -i /scratch2/10024/zhaochun/projects/farms_cdms_04192026/pulsepower/3d_lab_experiment/phase_field_hydromech/pulse_energy_tests/case_E5J/elasticity_mesh2x.i --allow-unused --split-mesh 500 --split-file /scratch2/10024/zhaochun/projects/farms_cdms_04192026/pulsepower/3d_lab_experiment/phase_field_hydromech/pulse_energy_tests/case_E5J/foo.cpr
ibrun /scratch2/10024/zhaochun/projects/farms_cdms_04192026/farms-opt -i /scratch2/10024/zhaochun/projects/farms_cdms_04192026/pulsepower/3d_lab_experiment/phase_field_hydromech/pulse_energy_tests/case_E5J/elasticity_mesh2x.i --allow-unused --use-split --split-file /scratch2/10024/zhaochun/projects/farms_cdms_04192026/pulsepower/3d_lab_experiment/phase_field_hydromech/pulse_energy_tests/case_E5J/foo.cpr --recover
