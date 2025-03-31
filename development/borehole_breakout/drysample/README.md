#!/bin/bash
#PBS -N drysample
#PBS -l select=5:ncpus=48:mpiprocs=20
#PBS -l walltime=72:00:00
#PBS -k doe
#PBS -P moose

cd $PBS_O_WORKDIR

module load use.moose moose-dev-openmpi/2025.03.03
mpiexec -n 100 moose-dev-exec ./farms-opt -i /scratch/zhaochun/projects/farms_cdbm_implicit/development/borehole_breakout/drysample/code_drysample_smallstrain.i --split-mesh 100 --split-file /scratch/zhaochun/projects/farms_cdbm_implicit/development/borehole_breakout/drysample/foo.cpr
mpiexec -n 100 moose-dev-exec ./farms-opt -i /scratch/zhaochun/projects/farms_cdbm_implicit/development/borehole_breakout/drysample/code_drysample_smallstrain.i --use-split --split-file /scratch/zhaochun/projects/farms_cdbm_implicit/development/borehole_breakout/drysample/foo.cpr