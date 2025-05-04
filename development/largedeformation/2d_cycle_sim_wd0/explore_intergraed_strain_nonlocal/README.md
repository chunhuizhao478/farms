Intergrated-Based Nonlocal Equvialent Strain Model

Segfault!

#!/bin/bash
#PBS -N integrated
#PBS -l select=4:ncpus=48:mpiprocs=30
#PBS -l walltime=72:00:00
#PBS -k doe
#PBS -P moose

cd $PBS_O_WORKDIR

module load use.moose moose-dev-openmpi/2025.04.22
mpiexec -n 120 moose-dev-exec ./farms-opt -i /scratch/zhaochun/projects/farms_cdbm_implicit/development/largedeformation/2d_cycle_sim_wd0/explore_intergraed_strain_nonlocal/static_solve_small/static_solve.i
mpiexec -n 120 moose-dev-exec ./farms-opt -i /scratch/zhaochun/projects/farms_cdbm_implicit/development/largedeformation/2d_cycle_sim_wd0/explore_intergraed_strain_nonlocal/dynamic_solve_small/dynamic_solve_main.i