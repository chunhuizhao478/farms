#!/bin/bash
#SBATCH -J cmame_2d_num_init_6fracture         # Job name
#SBATCH -o cmame_2d_num_init_6fracture.o%j     # Name of stdout output file
#SBATCH -e cmame_2d_num_init_6fracture.e%j     # Name of stderr error file
#SBATCH -p normal                              # Queue (partition) name
#SBATCH -N 4                                   # Total # of nodes
#SBATCH -n 200                                 # Total # of mpi tasks
#SBATCH -t 24:00:00                            # Run time (hh:mm:ss)
#SBATCH --mail-type=all                        # Send email at begin and end of job
#SBATCH -A EAR20006                            # Project/Allocation name
#SBATCH --mail-user=chunhui3@illinois.edu


# Load necessary modules
ml reset
ml gcc/11.2.0
ml impi/19.0.9
ml cuda/12.0
ml eigen/3.4.0
ml hdf5/1.14.6
ml netcdf/4.9.2
ml cmake/4.1.1

echo $CC $CXX $FC $F90 $F77
export CC=mpicc CXX=mpicxx FC=mpif90 F90=mpif90 F77=mpif77

export MOOSE_DIR=/work/10024/zhaochun/ls6/projects/moose-src
export PETSC_DIR=$MOOSE_DIR/petsc
export PETSC_ARCH=arch-moose

# Run the simulation
ibrun /scratch/10024/zhaochun/projects/farms_cdms/farms-opt -i /scratch/10024/zhaochun/projects/farms_cdms/pulsepower/cmame_revision/2d_solid/num_initial_fracture/case_6fracture/elasticity.i --allow-unused
