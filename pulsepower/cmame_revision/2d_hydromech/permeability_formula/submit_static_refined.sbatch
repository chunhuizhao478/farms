#!/bin/bash
#SBATCH -J cmame_2dhm_perm_static_ref          # Job name
#SBATCH -o cmame_2dhm_perm_static_ref.o%j      # Name of stdout output file
#SBATCH -e cmame_2dhm_perm_static_ref.e%j      # Name of stderr error file
#SBATCH -p normal                              # Queue (partition) name
#SBATCH -N 2                                   # Total # of nodes
#SBATCH -n 112                                 # Total # of mpi tasks
#SBATCH -t 04:00:00                            # Run time (hh:mm:ss)
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

# Refined-mesh static pre-solve. Produces static_solve_refined_out.e used as the
# initial condition for the *_refined.i elasticity runs (SolutionUserObject reads
# ../static_solve_refined_out.e). Submit this BEFORE any *_refined elasticity job.
# NOTE: 2d_mesh_refined.msh has lc = 5e-5 (vs 2e-4 baseline) => ~16x more elements;
# nodes/walltime bumped vs submit_static.sbatch. Increase further if it times out.
ibrun /scratch2/10024/zhaochun/projects/farms_cdms_04192026/farms-opt -i /scratch2/10024/zhaochun/projects/farms_cdms_04192026/pulsepower/cmame_revision/2d_hydromech/permeability_formula/static_solve_refined.i --allow-unused
