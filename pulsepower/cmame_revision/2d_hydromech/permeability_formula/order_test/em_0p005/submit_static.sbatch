#!/bin/bash
#SBATCH -J cmame_2dhm_ord2_0p005_static        # Job name
#SBATCH -o cmame_2dhm_ord2_0p005_static.o%j    # Name of stdout output file
#SBATCH -e cmame_2dhm_ord2_0p005_static.e%j    # Name of stderr error file
#SBATCH -p normal                              # Queue (partition) name
#SBATCH -N 1                                   # Total # of nodes
#SBATCH -n 56                                  # Total # of mpi tasks
#SBATCH -t 02:00:00                            # Run time (hh:mm:ss)
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

# Run the local static pre-solve (P1-P1 on 2d_mesh_o2.msh). Produces
# static_solve_out.e in this directory, which elasticity_E1d25.i reads via
# SolutionUserObject (mesh = ./static_solve_out.e). Using a local static on
# the same 2nd-order mesh avoids the 'point not found' errors that occur
# when interpolating the parent's P1-mesh static output at P2 DOFs of this
# case's 2nd-order mesh (the two mesh triangulations differ slightly near
# boundaries even though they come from the same .geo).
ibrun /scratch2/10024/zhaochun/projects/farms_cdms_04192026/farms-opt -i /scratch2/10024/zhaochun/projects/farms_cdms_04192026/pulsepower/cmame_revision/2d_hydromech/permeability_formula/order_test/em_0p005/static_solve.i --allow-unused
