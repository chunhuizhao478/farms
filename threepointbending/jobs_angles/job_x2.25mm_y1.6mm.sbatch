#!/bin/bash
#SBATCH -J 3dptrbend_whole_x2.25mm_y1.6mm        # Job name
#SBATCH -o 3dptrbend_whole_x2.25mm_y1.6mm.o%j    # Name of stdout output file
#SBATCH -e 3dptrbend_whole_x2.25mm_y1.6mm.e%j    # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 10               # Total # of nodes
#SBATCH -n 400              # Total # of mpi tasks
#SBATCH -t 48:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A EAR20006        # Project/Allocation name (req'd if you have more than 1)
#SBATCH --mail-user=chunhui3@illinois.edu

#=== USER PARAMETERS - MODIFY THIS ===
CASE_NAME="x2.25mm_y1.6mm"
#=====================================

# Derived paths
PROJECT_BASE="/scratch2/10024/zhaochun/projects/farms_cdms_01292026"
CASE_PATH="${PROJECT_BASE}/threepointbending/code_angles/case_whole_3d_${CASE_NAME}"
INPUT_PATH="${CASE_PATH}/elasticity.i"
SPLIT_PATH="${CASE_PATH}/foo.cpr"
NUM_SPLITS=400

# Load necessary modules
module swap intel gcc
module load cuda
export CXXFLAGS=-I/opt/apps/gcc/9.1.0/include/c++/9.1.0/
export CC=mpicc CXX=mpicxx FC=mpif90 F90=mpif90 F77=mpif77

# Enable MPI debugging
export MV2_DEBUG=1
export MV2_SHOW_ENV_INFO=1

export MOOSE_JOBS=6 METHODS=opt

# Split mesh
ibrun -n 1 ../../farms-opt -i ${INPUT_PATH} --split-mesh ${NUM_SPLITS} --split-file ${SPLIT_PATH}

# Run simulation
ibrun ../../farms-opt -i ${INPUT_PATH} --use-split --split-file ${SPLIT_PATH} --allow-unused
