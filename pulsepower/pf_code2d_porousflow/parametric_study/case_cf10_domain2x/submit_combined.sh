#!/bin/bash
#SBATCH -J case_cf10_domain2x_combined        # Job name
#SBATCH -o case_cf10_domain2x_combined.o%j    # Name of stdout output file
#SBATCH -e case_cf10_domain2x_combined.e%j    # Name of stderr error file
#SBATCH -p normal       # Queue (partition) name
#SBATCH -N 8           # Total # of nodes
#SBATCH -n 200          # Total # of mpi tasks
#SBATCH -t 48:00:00            # Run time (hh:mm:ss)
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A EAR20006         # Project/Allocation name (req'd if you have more than 1)
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

# Change to case directory
cd /scratch1/10024/zhaochun/projects/farms_cdms_11022025/pulsepower/pf_code2d_porousflow/parametric_study/case_cf10_domain2x

echo "================================================================================"
echo "Starting combined static + dynamic simulation"
echo "Case: case_cf10_domain2x"
echo "Start time: $(date)"
echo "================================================================================"

# Step 1: Run static solve
echo ""
echo "--------------------------------------------------------------------------------"
echo "STEP 1: Running static solve..."
echo "--------------------------------------------------------------------------------"
echo "Input file: /scratch1/10024/zhaochun/projects/farms_cdms_11022025/pulsepower/pf_code2d_porousflow/parametric_study/case_cf10_domain2x/static_solve.i"
echo "Start time: $(date)"
echo ""

ibrun ./farms-opt -i /scratch1/10024/zhaochun/projects/farms_cdms_11022025/pulsepower/pf_code2d_porousflow/parametric_study/case_cf10_domain2x/static_solve.i --allow-unused

STATIC_EXIT_CODE=$?

echo ""
echo "Static solve completed with exit code: $STATIC_EXIT_CODE"
echo "End time: $(date)"

if [ $STATIC_EXIT_CODE -ne 0 ]; then
    echo ""
    echo "ERROR: Static solve failed with exit code $STATIC_EXIT_CODE"
    echo "Aborting job..."
    exit $STATIC_EXIT_CODE
fi

# Check if static solve output exists
if [ ! -f "static_solve_out.e" ]; then
    echo ""
    echo "ERROR: Static solve output file (static_solve_out.e) not found!"
    echo "Aborting job..."
    exit 1
fi

echo ""
echo "Static solve output verified: static_solve_out.e exists"

# Step 2: Run dynamic solve (elasticity)
echo ""
echo "--------------------------------------------------------------------------------"
echo "STEP 2: Running dynamic solve (elasticity)..."
echo "--------------------------------------------------------------------------------"
echo "Input file: /scratch1/10024/zhaochun/projects/farms_cdms_11022025/pulsepower/pf_code2d_porousflow/parametric_study/case_cf10_domain2x/elasticity.i"
echo "Start time: $(date)"
echo ""

ibrun ./farms-opt -i /scratch1/10024/zhaochun/projects/farms_cdms_11022025/pulsepower/pf_code2d_porousflow/parametric_study/case_cf10_domain2x/elasticity.i --allow-unused

DYNAMIC_EXIT_CODE=$?

echo ""
echo "Dynamic solve completed with exit code: $DYNAMIC_EXIT_CODE"
echo "End time: $(date)"

if [ $DYNAMIC_EXIT_CODE -ne 0 ]; then
    echo ""
    echo "ERROR: Dynamic solve failed with exit code $DYNAMIC_EXIT_CODE"
    exit $DYNAMIC_EXIT_CODE
fi

echo ""
echo "================================================================================"
echo "Combined simulation completed successfully!"
echo "End time: $(date)"
echo "================================================================================"

exit 0
