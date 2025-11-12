#!/bin/bash
#SBATCH -J case_cf5_domain5x_elasticity        # Job name
#SBATCH -o case_cf5_domain5x_elasticity.o%j    # Name of stdout output file
#SBATCH -e case_cf5_domain5x_elasticity.e%j    # Name of stderr error file
#SBATCH -p normal       # Queue (partition) name
#SBATCH -N 4           # Total # of nodes
#SBATCH -n 200          # Total # of mpi tasks
#SBATCH -t 24:00:00            # Run time (hh:mm:ss)
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A ASC25056         # Project/Allocation name (req'd if you have more than 1)
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
ibrun ./farms-opt -i /scratch/10024/zhaochun/projects/farms_cdms/pulsepower/pf_code2d_puresolid/parametric_study/case_cf5_domain5x/elasticity.i --allow-unused
