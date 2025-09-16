#!/bin/bash
#SBATCH -J tpv2053d        # Job name
#SBATCH -o tpv2053d.o%j    # Name of stdout output file
#SBATCH -e tpv2053d.e%j    # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 6               # Total # of nodes 
#SBATCH -n 100             # Total # of mpi tasks
#SBATCH -t 48:00:00        # Run time (hh:mm:ss)
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

ibrun ./farms-opt -i /scratch1/10024/zhaochun/projects/farms_benchmark/examples/benchmark_tpv2053D/czm/tpv2053D_tet.i