#!/usr/bin/env python3
"""
Generate Combined SLURM Submission Scripts

This script generates a single sbatch file per case that runs:
1. Static solve (static_solve.i) first
2. Dynamic solve (elasticity.i) second

Usage:
    python3 generate_combined_submit.py

Author: Auto-generated
Date: 2025-11-09
"""

import sys
from pathlib import Path

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Base directory (where this script is located)
BASE_DIR = Path(__file__).parent

# Project directory on HPC (ABSOLUTE PATH on the cluster)
# This should be the path to your project root on the HPC system
HPC_PROJECT_ROOT = "/scratch/10024/zhaochun/projects/farms_cdms"

# Relative path from project root to this parametric study directory
RELATIVE_STUDY_PATH = "pulsepower/pf_code2d_porousflow/parametric_study"

# Executable name
EXECUTABLE = "/scratch/10024/zhaochun/projects/farms_cdms/farms-opt"

# Input file to run (elasticity.i or static_solve.i)
INPUT_FILE = "elasticity.i"

# SLURM Parameters (defaults)
SLURM_PARAMS = {
    "partition": "normal",  # Queue name
    "nodes": 4,  # Number of nodes
    "ntasks": 200,  # Total number of MPI tasks
    "time": "24:00:00",  # Wall time (hh:mm:ss)
    "account": "ASC25096",  # Project/Allocation name
    "mail_user": "chunhui3@illinois.edu",  # Email for notifications
    "mail_type": "all",  # Email notification type
}

# Module loading commands (customize for your HPC system)
MODULE_COMMANDS = """# Load necessary modules
ml reset
ml gcc/11.2.0
ml impi/19.0.9
ml cuda/12.0
ml eigen/3.4.0
ml hdf5/1.14.6
ml netcdf/4.9.2
ml cmake/4.1.1

echo $CC $CXX $FC $F90 $F77
export CC=mpicc CXX=mpicxx FC=mpif90 F90=mpif90 F77=mpif77

export MOOSE_DIR=/work/10024/zhaochun/ls6/projects/moose-src
export PETSC_DIR=$MOOSE_DIR/petsc
export PETSC_ARCH=arch-moose"""

# ==============================================================================
# SLURM JOB TEMPLATE
# ==============================================================================

COMBINED_TEMPLATE = """#!/bin/bash
#SBATCH -J {job_name}        # Job name
#SBATCH -o {job_name}.o%j    # Name of stdout output file
#SBATCH -e {job_name}.e%j    # Name of stderr error file
#SBATCH -p {partition}       # Queue (partition) name
#SBATCH -N {nodes}           # Total # of nodes
#SBATCH -n {ntasks}          # Total # of mpi tasks
#SBATCH -t {time}            # Run time (hh:mm:ss)
#SBATCH --mail-type={mail_type}    # Send email at begin and end of job
#SBATCH -A {account}         # Project/Allocation name (req'd if you have more than 1)
#SBATCH --mail-user={mail_user}

{module_commands}

# Change to case directory
cd {case_dir}

echo "================================================================================"
echo "Starting combined static + dynamic simulation"
echo "Case: {case_name}"
echo "Start time: $(date)"
echo "================================================================================"

# Step 1: Run static solve
echo ""
echo "--------------------------------------------------------------------------------"
echo "STEP 1: Running static solve..."
echo "--------------------------------------------------------------------------------"
echo "Input file: {static_input}"
echo "Start time: $(date)"
echo ""

ibrun {executable} -i {static_input} --allow-unused

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
echo "Input file: {dynamic_input}"
echo "Start time: $(date)"
echo ""

ibrun {executable} -i {dynamic_input} --allow-unused

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
"""

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================


def find_case_folders():
    """Find all case folders in the parametric study directory."""
    case_folders = sorted(BASE_DIR.glob("case_*"))
    case_folders = [f for f in case_folders if f.is_dir()]
    return case_folders


def case_has_static_solve(case_folder):
    """Check if a case has static_solve.i file."""
    return (case_folder / "static_solve.i").exists()


def generate_combined_script(case_folder):
    """Generate combined SLURM script content for a case."""
    case_name = case_folder.name

    # Construct absolute paths on HPC
    case_dir = f"{HPC_PROJECT_ROOT}/{RELATIVE_STUDY_PATH}/{case_name}"
    static_input = f"{case_dir}/static_solve.i"
    dynamic_input = f"{case_dir}/elasticity.i"

    # Job name
    job_name = f"{case_name}_combined"

    # Fill in the template
    script_content = COMBINED_TEMPLATE.format(
        job_name=job_name,
        partition=SLURM_PARAMS["partition"],
        nodes=SLURM_PARAMS["nodes"],
        ntasks=SLURM_PARAMS["ntasks"],
        time=SLURM_PARAMS["time"],
        mail_type=SLURM_PARAMS["mail_type"],
        account=SLURM_PARAMS["account"],
        mail_user=SLURM_PARAMS["mail_user"],
        module_commands=MODULE_COMMANDS,
        case_dir=case_dir,
        case_name=case_name,
        static_input=static_input,
        dynamic_input=dynamic_input,
        executable=EXECUTABLE,
    )

    return script_content


def write_combined_script(case_folder, script_content):
    """Write combined SLURM script to file."""
    script_filename = "submit_combined.sh"
    script_path = case_folder / script_filename

    # Write the script
    with open(script_path, "w") as f:
        f.write(script_content)

    # Make it executable
    script_path.chmod(0o755)

    return script_path


# ==============================================================================
# MAIN EXECUTION
# ==============================================================================


def main():
    print("=" * 80)
    print("Combined SLURM Job Script Generator")
    print("=" * 80)
    print(f"\nBase directory: {BASE_DIR}")
    print(f"HPC project root: {HPC_PROJECT_ROOT}")
    print(f"Study path: {RELATIVE_STUDY_PATH}")

    # Find all case folders
    case_folders = find_case_folders()

    if not case_folders:
        print("\n✗ ERROR: No case folders found")
        sys.exit(1)

    print(f"\nFound {len(case_folders)} case(s)")

    # Separate cases by whether they have static solve
    cases_with_static = []
    cases_without_static = []

    for case_folder in case_folders:
        if case_has_static_solve(case_folder):
            cases_with_static.append(case_folder)
        else:
            cases_without_static.append(case_folder)

    print(f"\nCases with static solve: {len(cases_with_static)}")
    for cf in cases_with_static:
        print(f"  - {cf.name}")

    if cases_without_static:
        print(
            f"\nCases WITHOUT static solve (will be skipped): {len(cases_without_static)}"
        )
        for cf in cases_without_static:
            print(f"  - {cf.name}")

    if not cases_with_static:
        print("\n✗ ERROR: No cases with static_solve.i found")
        sys.exit(1)

    # Generate scripts
    print("\n" + "=" * 80)
    print("Generating combined SLURM job scripts...")
    print("=" * 80)

    generated_scripts = []
    for case_folder in cases_with_static:
        # Generate script content
        script_content = generate_combined_script(case_folder)

        # Write script to file
        script_path = write_combined_script(case_folder, script_content)

        generated_scripts.append(script_path)
        print(f"  ✓ Created: {script_path.relative_to(BASE_DIR)}")

    # Print summary
    print("\n" + "=" * 80)
    print(f"✓ Successfully generated {len(generated_scripts)} combined SLURM scripts!")
    print("=" * 80)

    print("\nGenerated scripts:")
    for script_path in generated_scripts:
        print(f"  - {script_path.relative_to(BASE_DIR)}")

    # Print usage instructions
    print("\n" + "=" * 80)
    print("Usage Instructions:")
    print("=" * 80)

    print("\n1. Transfer files to HPC cluster")
    print(f"   (files should already be on cluster)")

    print("\n2. Submit a single job:")
    print(f"   cd {HPC_PROJECT_ROOT}/{RELATIVE_STUDY_PATH}/{cases_with_static[0].name}")
    print("   sbatch submit_combined.sh")

    print("\n3. Submit all jobs:")
    print(f"   cd {HPC_PROJECT_ROOT}/{RELATIVE_STUDY_PATH}")
    print("   for dir in case_*/; do")
    print('     if [ -f "$dir/submit_combined.sh" ]; then')
    print('       cd "$dir" && sbatch submit_combined.sh && cd ..')
    print("     fi")
    print("   done")

    print("\n4. Monitor jobs:")
    print("   squeue -u $USER")

    print("\n5. Check output:")
    print("   tail -f case_*/*_combined.o*")

    print("\n" + "=" * 80)
    print("Notes:")
    print("  - Each script runs static_solve.i THEN elasticity.i sequentially")
    print("  - The script will abort if static solve fails")
    print("  - Both steps run in the same SLURM allocation")
    print("  - Total wall time is the same as individual jobs (48:00:00)")
    print("=" * 80)


if __name__ == "__main__":
    main()
