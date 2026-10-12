#!/bin/bash
#SBATCH -J 3dptrbend_dtprobe_p1        # Job name
#SBATCH -o 3dptrbend_dtprobe_p1.o%j    # Name of stdout output file
#SBATCH -e 3dptrbend_dtprobe_p1.e%j    # Name of stderr error file
#SBATCH -p development     # Queue (partition) name
#SBATCH -N 20               # Total # of nodes
#SBATCH -n 800              # Total # of mpi tasks
#SBATCH -t 02:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A EAR20006        # Project/Allocation name (req'd if you have more than 1)
#SBATCH --mail-user=chunhui3@illinois.edu

# dt probes (dt 0.1 to t = 50, make_dt_probe.py) of 2 cases side by side, each on its own
# 10 nodes / 400 ranks (ibrun -n 400 -o offset). Each probe writes split.log,
# probe.log, elasticity_out.csv and elasticity_exodus.e in its own case folder.
#=== USER PARAMETERS - MODIFY THIS ===
CASES=(
  "code_2holes_connected/case_2holes_connected_3d_DL2.0mm_DR2.0mm_x2.00mm_y3.6mm_upline_dt0p1_end50"
  "code_pillars/case_pillars_3d_welded_mismatch_2onN_dt0p1_end50"
)
#=====================================

# Derived paths
PROJECT_BASE="/scratch2/10024/zhaochun/projects/farms_cdms_04192026"
RANKS_PER_CASE=400

# Load necessary modules
module swap intel gcc
module load cuda
export CXXFLAGS=-I/opt/apps/gcc/9.1.0/include/c++/9.1.0/
export CC=mpicc CXX=mpicxx FC=mpif90 F90=mpif90 F77=mpif77

# Enable MPI debugging
export MV2_DEBUG=1
export MV2_SHOW_ENV_INFO=1

export MOOSE_JOBS=6 METHODS=opt

# Split each mesh (serial), side by side
for i in "${!CASES[@]}"; do
  C="${PROJECT_BASE}/threepointbending/${CASES[$i]}"
  ibrun -n 1 -o $((i * RANKS_PER_CASE)) ../../farms-opt -i ${C}/elasticity.i --split-mesh ${RANKS_PER_CASE} --split-file ${C}/foo.cpr > ${C}/split.log 2>&1 &
done
wait

# Run each probe on its own nodes, side by side
for i in "${!CASES[@]}"; do
  C="${PROJECT_BASE}/threepointbending/${CASES[$i]}"
  ibrun -n ${RANKS_PER_CASE} -o $((i * RANKS_PER_CASE)) ../../farms-opt -i ${C}/elasticity.i --use-split --split-file ${C}/foo.cpr --allow-unused > ${C}/probe.log 2>&1 &
done
wait
