#!/usr/bin/env python3
"""
SLURM Job Generator for Porousflow Parametric Studies

This script generates SLURM job submission files for all parametric study cases.
Supports both static_solve.i and elasticity.i simulations.

Usage:
    # Generate elasticity.i job files for all cases
    python3 generate_slurm_jobs.py

    # Generate static_solve.i job files
    python3 generate_slurm_jobs.py --input-file static_solve.i

    # Generate for specific pattern
    python3 generate_slurm_jobs.py --pattern "case_cf*"

    # Custom configuration
    python3 generate_slurm_jobs.py --nodes 4 --ntasks 100 --time 24:00:00

Author: Auto-generated
Date: 2025-11-08
"""

import argparse
import sys
from pathlib import Path

# ==============================================================================
# CONFIGURATION - Modify these parameters as needed
# ==============================================================================

# Base directory (where this script is located)
BASE_DIR = Path(__file__).parent

# Project directory on HPC (ABSOLUTE PATH on the cluster)
# This should be the path to your project root on the HPC system
HPC_PROJECT_ROOT = "/scratch1/10024/zhaochun/projects/farms_cdms_11022025"

# Relative path from project root to this parametric study directory
RELATIVE_STUDY_PATH = "pulsepower/pf_code2d_porousflow/parametric_study"

# Executable name
EXECUTABLE = "./farms-opt"

# Input file to run (elasticity.i or static_solve.i)
INPUT_FILE = "elasticity.i"

# SLURM Parameters (defaults)
SLURM_DEFAULTS = {
    "partition": "normal",  # Queue name
    "nodes": 4,  # Number of nodes
    "ntasks": 200,  # Total number of MPI tasks
    "time": "24:00:00",  # Wall time (hh:mm:ss)
    "account": "ASC25056",  # Project/Allocation name
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

# Additional SLURM flags (optional)
ADDITIONAL_FLAGS = [
    # Add any additional SLURM flags here, e.g.:
    # '--constraint=haswell',
    # '--exclusive',
]

# ==============================================================================
# SLURM JOB TEMPLATE
# ==============================================================================

SLURM_TEMPLATE = """#!/bin/bash
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
{additional_flags}

{module_commands}

# Run the simulation
ibrun {executable} -i {input_path} --allow-unused
"""

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================


def find_case_folders(pattern=None):
    """
    Find case folders in the parametric study directory.

    Args:
        pattern (str): Glob pattern to match case folders (e.g., "case_cf*")
                      If None, returns all case_* folders

    Returns:
        list: List of Path objects for case folders
    """
    if pattern is None:
        pattern = "case_*"

    case_folders = sorted(BASE_DIR.glob(pattern))
    # Filter to only directories
    case_folders = [f for f in case_folders if f.is_dir()]

    return case_folders


def generate_job_name(case_name, input_file):
    """
    Generate a descriptive job name from case folder name.

    Args:
        case_name (str): Name of the case folder
        input_file (str): Input file being run

    Returns:
        str: Job name
    """
    # Remove .i extension from input file
    input_prefix = input_file.replace(".i", "")

    # Create job name
    job_name = f"{case_name}_{input_prefix}"

    return job_name


def generate_slurm_script(case_folder, input_file, slurm_params):
    """
    Generate SLURM job script content for a case.

    Args:
        case_folder (Path): Path to the case folder
        input_file (str): Input file to run (e.g., 'elasticity.i')
        slurm_params (dict): SLURM parameters

    Returns:
        str: SLURM script content
    """
    case_name = case_folder.name

    # Generate job name
    job_name = generate_job_name(case_name, input_file)

    # Construct absolute path to input file on HPC
    input_path = f"{HPC_PROJECT_ROOT}/{RELATIVE_STUDY_PATH}/{case_name}/{input_file}"

    # Format additional flags
    additional_flags = ""
    if ADDITIONAL_FLAGS:
        additional_flags = "\n".join(f"#SBATCH {flag}" for flag in ADDITIONAL_FLAGS)

    # Fill in the template
    script_content = SLURM_TEMPLATE.format(
        job_name=job_name,
        partition=slurm_params["partition"],
        nodes=slurm_params["nodes"],
        ntasks=slurm_params["ntasks"],
        time=slurm_params["time"],
        mail_type=slurm_params["mail_type"],
        account=slurm_params["account"],
        mail_user=slurm_params["mail_user"],
        additional_flags=additional_flags,
        module_commands=MODULE_COMMANDS,
        executable=EXECUTABLE,
        input_path=input_path,
    )

    return script_content


def write_slurm_script(case_folder, script_content, input_file):
    """
    Write SLURM script to file.

    Args:
        case_folder (Path): Path to the case folder
        script_content (str): SLURM script content
        input_file (str): Input file being run

    Returns:
        Path: Path to the written script file
    """
    # Generate script filename
    input_prefix = input_file.replace(".i", "")
    script_filename = f"submit_{input_prefix}.sh"
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
    parser = argparse.ArgumentParser(
        description="Generate SLURM job submission scripts for parametric cases"
    )

    # Case selection
    case_group = parser.add_mutually_exclusive_group()
    case_group.add_argument(
        "--case",
        type=str,
        help="Generate for a specific case folder (e.g., case_cf1_domain1x)",
    )
    case_group.add_argument(
        "--pattern",
        type=str,
        help='Generate for all cases matching pattern (e.g., "case_cf*")',
    )

    # SLURM parameters
    parser.add_argument(
        "--nodes",
        type=int,
        default=SLURM_DEFAULTS["nodes"],
        help=f"Number of nodes (default: {SLURM_DEFAULTS['nodes']})",
    )
    parser.add_argument(
        "--ntasks",
        "-n",
        type=int,
        default=SLURM_DEFAULTS["ntasks"],
        help=f"Total number of MPI tasks (default: {SLURM_DEFAULTS['ntasks']})",
    )
    parser.add_argument(
        "--time",
        "-t",
        type=str,
        default=SLURM_DEFAULTS["time"],
        help=f"Wall time in hh:mm:ss (default: {SLURM_DEFAULTS['time']})",
    )
    parser.add_argument(
        "--partition",
        "-p",
        type=str,
        default=SLURM_DEFAULTS["partition"],
        help=f"Queue/partition name (default: {SLURM_DEFAULTS['partition']})",
    )
    parser.add_argument(
        "--account",
        "-A",
        type=str,
        default=SLURM_DEFAULTS["account"],
        help=f"Project/Allocation name (default: {SLURM_DEFAULTS['account']})",
    )
    parser.add_argument(
        "--input-file",
        "-i",
        type=str,
        default=INPUT_FILE,
        help=f"Input file to run (default: {INPUT_FILE})",
    )

    args = parser.parse_args()

    # Update SLURM parameters from command line
    slurm_params = SLURM_DEFAULTS.copy()
    slurm_params["nodes"] = args.nodes
    slurm_params["ntasks"] = args.ntasks
    slurm_params["time"] = args.time
    slurm_params["partition"] = args.partition
    slurm_params["account"] = args.account

    # Print header
    print("=" * 80)
    print("SLURM Job Script Generator - Porousflow Parametric Study")
    print("=" * 80)
    print(f"\nBase directory: {BASE_DIR}")
    print(f"HPC project root: {HPC_PROJECT_ROOT}")
    print(f"Study path: {RELATIVE_STUDY_PATH}")
    print(f"Input file: {args.input_file}")

    # Determine which cases to process
    if args.case:
        case_folder = BASE_DIR / args.case
        if not case_folder.is_dir():
            print(f"\n✗ ERROR: Case folder not found: {case_folder}")
            sys.exit(1)
        case_folders = [case_folder]
    elif args.pattern:
        case_folders = find_case_folders(args.pattern)
        if not case_folders:
            print(f"\n✗ ERROR: No case folders found matching pattern: {args.pattern}")
            sys.exit(1)
    else:  # All cases
        case_folders = find_case_folders()
        if not case_folders:
            print(f"\n✗ ERROR: No case folders found")
            sys.exit(1)

    print(f"\nFound {len(case_folders)} case(s) to process:")
    for cf in case_folders:
        print(f"  - {cf.name}")

    print("\nSLURM Parameters:")
    print(f"  Nodes: {slurm_params['nodes']}")
    print(f"  Tasks: {slurm_params['ntasks']}")
    print(f"  Time: {slurm_params['time']}")
    print(f"  Partition: {slurm_params['partition']}")
    print(f"  Account: {slurm_params['account']}")

    # Generate scripts
    print("\n" + "=" * 80)
    print("Generating SLURM job scripts...")
    print("=" * 80)

    generated_scripts = []
    for case_folder in case_folders:
        # Generate script content
        script_content = generate_slurm_script(
            case_folder, args.input_file, slurm_params
        )

        # Write script to file
        script_path = write_slurm_script(case_folder, script_content, args.input_file)

        generated_scripts.append(script_path)
        print(f"  ✓ Created: {script_path.relative_to(BASE_DIR)}")

    # Print summary
    print("\n" + "=" * 80)
    print(f"✓ Successfully generated {len(generated_scripts)} SLURM job scripts!")
    print("=" * 80)

    print("\nGenerated scripts:")
    for script_path in generated_scripts:
        print(f"  - {script_path.relative_to(BASE_DIR)}")

    # Print submission instructions
    print("\n" + "=" * 80)
    print("Workflow for Porousflow Simulations:")
    print("=" * 80)

    if args.input_file == "static_solve.i":
        print("\n📋 STEP 1: Run static_solve.i to get initial conditions")
        print("   1a. Transfer files to HPC")
        print("   1b. Submit static solve jobs:")
        print(f"       cd {HPC_PROJECT_ROOT}/{RELATIVE_STUDY_PATH}")
        print(
            f'       for dir in case_*/; do cd "$dir" && sbatch submit_static_solve.sh && cd ..; done'
        )
        print("   1c. Wait for completion and verify results")
        print("   1d. Run update script locally to extract energy values:")
        print("       python3 run_static_solve_and_update.py --all --skip-simulation")
    else:
        print("\n📋 STEP 2: Run elasticity.i for dynamic simulation")
        print("   2a. Ensure static_solve.i has completed and elasticity.i is updated")
        print("   2b. Transfer updated files to HPC")
        print("   2c. Submit elasticity jobs:")
        print(f"       cd {HPC_PROJECT_ROOT}/{RELATIVE_STUDY_PATH}")
        print(
            f'       for dir in case_*/; do cd "$dir" && sbatch submit_elasticity.sh && cd ..; done'
        )

    print("\n" + "=" * 80)
    print("Submission Commands:")
    print("=" * 80)
    print("\n1. Transfer files to HPC:")
    print(
        f"   scp -r {BASE_DIR} username@frontera.tacc.utexas.edu:{HPC_PROJECT_ROOT}/{RELATIVE_STUDY_PATH}"
    )

    print("\n2. Submit jobs on HPC:")
    print("   # Submit a single job")
    print(f"   cd {HPC_PROJECT_ROOT}/{RELATIVE_STUDY_PATH}/case_cf1_domain1x")
    print(f"   sbatch submit_{args.input_file.replace('.i', '')}.sh")

    print("\n   # Submit all jobs")
    print(f"   cd {HPC_PROJECT_ROOT}/{RELATIVE_STUDY_PATH}")
    print(
        f'   for dir in case_*/; do cd "$dir" && sbatch submit_{args.input_file.replace(".i", "")}.sh && cd ..; done'
    )

    print("\n3. Monitor jobs:")
    print("   squeue -u $USER")

    print("\n4. Check job output:")
    print("   # Look for .o and .e files in each case directory")
    print("   tail -f case_cf1_domain1x/*.o*")

    print("\n" + "=" * 80)
    print("Notes:")
    print(f"  - Edit CONFIGURATION section in this script to customize parameters")
    print(f"  - Update HPC_PROJECT_ROOT to match your actual path on the cluster")
    print(f"  - For porousflow: run static_solve.i first, then elasticity.i")
    print(f"  - Verify module commands match your HPC environment")
    print("=" * 80)


if __name__ == "__main__":
    main()
