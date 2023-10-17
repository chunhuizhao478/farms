#!/bin/bash
#PBS -N contmf
#PBS -l select=5:ncpus=48:mpiprocs=36
#PBS -l walltime=72:00:00
#PBS -P moose

cd $PBS_O_WORKDIR
source /etc/profile.d/modules.sh

module load use.moose moose-dev
mpiexec ./farms-opt -i /home/zhaochun/projects_jul2/farms/examples/cdbm_planarfault/contmf_main.i 
