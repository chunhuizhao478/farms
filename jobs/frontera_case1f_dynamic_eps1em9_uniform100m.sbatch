#!/bin/bash
#SBATCH -J case1f_dyn_eps1em9_u100m
#SBATCH -o case1f_dyn_eps1em9_u100m.o%j
#SBATCH -e case1f_dyn_eps1em9_u100m.e%j
#SBATCH -p normal
#SBATCH -N 12
#SBATCH -n 672
#SBATCH -t 48:00:00
#SBATCH --mail-type=all
#SBATCH -A EAR20006
#SBATCH --mail-user=chunhui3@illinois.edu

set -euo pipefail

module swap intel gcc
module swap impi mvapich2-x
module load cuda

export OMP_NUM_THREADS=1
export MOOSE_JOBS=6
export METHODS=opt

export CC=mpicc
export CXX=mpicxx
export FC=mpif90
export F90=mpif90
export F77=mpif77
export CXXFLAGS=-I/opt/apps/gcc/9.1.0/include/c++/9.1.0/

export MV2_DEBUG=1
export MV2_SHOW_ENV_INFO=1

APP_ROOT=/scratch2/10024/zhaochun/projects/farms_benchmark_03302026
INPUT_FILE=$APP_ROOT/development/paper_draft_depth_dependent/case1f/dynamic_solve/dynamic_solve_eps1em9_uniform100m.i

cd "$APP_ROOT"

# Step 1: Pre-split mesh on single rank (avoids full-mesh-on-every-rank OOM)
echo "=== Splitting mesh into 672 partitions ==="
ibrun -np 1 ./farms-opt --split-mesh 672 -i "$INPUT_FILE"

# Step 2: Run solve using pre-split mesh
echo "=== Running dynamic solve with pre-split mesh ==="
ibrun ./farms-opt -i "$INPUT_FILE" --use-split
