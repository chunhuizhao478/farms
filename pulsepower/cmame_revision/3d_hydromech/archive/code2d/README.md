#!/bin/bash
#SBATCH -J test1_2d_pf        # Job name
#SBATCH -o test1_2d_pf.o%j    # Name of stdout output file
#SBATCH -e test1_2d_pf.e%j    # Name of stderr error file
#SBATCH -p development           # Queue (partition) name
#SBATCH -N 20               # Total # of nodes 
#SBATCH -n 200              # Total # of mpi tasks
#SBATCH -t 02:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A EAR20006        # Project/Allocation name (req'd if you have more than 1)
#SBATCH --mail-user=chunhui3@illinois.edu

# Load necessary modules
module swap intel gcc
#module swap impi mvapich2-x
module load cuda
export CXXFLAGS=-I/opt/apps/gcc/9.1.0/include/c++/9.1.0/
export CC=mpicc CXX=mpicxx FC=mpif90 F90=mpif90 F77=mpif77

# Set compilers and flags
export CC=mpicc
export CXX=mpicxx
export FC=mpif90
export F90=mpif90
export F77=mpif77

export CXXFLAGS=-I/opt/apps/gcc/9.1.0/include/c++/9.1.0/
export CC=mpicc CXX=mpicxx FC=mpif90 F90=mpif90 F77=mpif77

# Enable MPI debugging
export MV2_DEBUG=1
export MV2_SHOW_ENV_INFO=1

export MOOSE_JOBS=6 METHODS=opt

# Define the list as an array

ibrun ./farms-opt -i /scratch1/10024/zhaochun/projects/farms_cdms/pulsepower/porousflowcoupling/code2d/static_solve.i
ibrun -n 1 ./farms-opt -i /scratch1/10024/zhaochun/projects/farms_cdms/pulsepower/porousflowcoupling/code2d/elasticity.i --split-mesh 200 --split-file /scratch1/10024/zhaochun/projects/farms_cdms/pulsepower/porousflowcoupling/code2d/foo.cpr
ibrun ./farms-opt -i /scratch1/10024/zhaochun/projects/farms_cdms/pulsepower/porousflowcoupling/code2d/elasticity.i --use-split --split-file /scratch1/10024/zhaochun/projects/farms_cdms/pulsepower/porousflowcoupling/code2d/foo.cpr