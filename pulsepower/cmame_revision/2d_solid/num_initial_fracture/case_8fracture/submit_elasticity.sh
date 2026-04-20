#!/bin/bash
#SBATCH -J cmame_2d_num_init_8fracture         # Job name
#SBATCH -o cmame_2d_num_init_8fracture.o%j     # Name of stdout output file
#SBATCH -e cmame_2d_num_init_8fracture.e%j     # Name of stderr error file
#SBATCH -p normal                              # Queue (partition) name
#SBATCH -N 4                                   # Total # of nodes
#SBATCH -n 200                                 # Total # of mpi tasks
#SBATCH -t 24:00:00                            # Run time (hh:mm:ss)
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
ibrun /scratch/10024/zhaochun/projects/farms_cdms/farms-opt -i /scratch/10024/zhaochun/projects/farms_cdms/pulsepower/cmame_revision/2d_solid/num_initial_fracture/case_8fracture/elasticity.i --allow-unused
