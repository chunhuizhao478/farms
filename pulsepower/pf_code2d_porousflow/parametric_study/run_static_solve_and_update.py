#!/usr/bin/env python3
"""
Static Solve Runner and Elasticity File Updater

This script:
1. Runs static_solve.i for specified case(s)
2. Reads the output CSV file (static_solve_out_csv.csv)
3. Extracts energy values from the CSV
4. Updates the elasticity.i file with the extracted values

Usage:
    # Run for a specific case
    python3 run_static_solve_and_update.py --case case_cf1_domain1x

    # Run for all cases matching a pattern
    python3 run_static_solve_and_update.py --pattern "case_cf*"

    # Run for all cases
    python3 run_static_solve_and_update.py --all

    # Dry run (extract values without running simulation)
    python3 run_static_solve_and_update.py --case case_cf1_domain1x --dry-run

Author: Auto-generated
Date: 2025-11-08
"""

import os
import sys
import argparse
import subprocess
import csv
from pathlib import Path
import glob
import time

# ==============================================================================
# CONFIGURATION - Modify these parameters as needed
# ==============================================================================

# Number of MPI processes to use
MPI_PROCESSES = 8

# Executable name (relative to project root or absolute path)
EXECUTABLE = "./farms-opt"

# Base directory (where this script is located)
SCRIPT_DIR = Path(__file__).parent

# Project root directory (where farms-opt executable is located)
# Assuming project structure: farms_cdms/pulsepower/pf_code2d_porousflow/parametric_study/
PROJECT_ROOT = SCRIPT_DIR.parent.parent.parent

# CSV output filename (default MOOSE output pattern)
CSV_FILENAME = "static_solve_csv.csv"

# Variables to extract from CSV (column names in CSV file)
CSV_COLUMN_NAMES = {
    'fluid_elastic_energy_total_static': 'fluid_elastic_energy_total_static',
    'solid_elastic_energy_total_static': 'solid_elastic_energy_static',  # Note: different name in CSV
    'full_input_energy_static': 'full_input_energy_static'
}

# Variables to write to elasticity.i (in order)
ELASTICITY_VARIABLES = [
    'fluid_elastic_energy_total_static',
    'solid_elastic_energy_total_static',
    'full_input_energy_static'
]

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

    case_folders = sorted(SCRIPT_DIR.glob(pattern))
    # Filter to only directories
    case_folders = [f for f in case_folders if f.is_dir()]

    return case_folders


def run_static_solve(case_folder, mpi_processes=8, quiet=False):
    """
    Run static_solve.i for a given case folder.

    Args:
        case_folder (Path): Path to the case folder
        mpi_processes (int): Number of MPI processes to use
        quiet (bool): If True, suppress MOOSE output (faster)

    Returns:
        bool: True if successful, False otherwise
    """
    static_solve_file = case_folder / "static_solve.i"

    if not static_solve_file.exists():
        print(f"  ✗ ERROR: static_solve.i not found in {case_folder.name}")
        return False

    # Construct the relative path from project root
    relative_input_path = static_solve_file.relative_to(PROJECT_ROOT)

    # Construct the mpirun command
    cmd = [
        "mpirun",
        "-np", str(mpi_processes),
        EXECUTABLE,
        "-i", str(relative_input_path)
    ]

    print(f"\n  Running: {' '.join(cmd)}")
    print(f"  Working directory: {PROJECT_ROOT}")

    if quiet:
        print(f"  Running in QUIET mode (output suppressed for speed)...")
    else:
        print(f"  MOOSE output will be displayed below...")
        print(f"  " + "-" * 76)

    try:
        # Run the command from the project root directory
        # For performance: let output go directly to terminal (or /dev/null if quiet)
        # This avoids buffering massive amounts of MOOSE output in memory
        if quiet:
            # Suppress output for speed
            result = subprocess.run(
                cmd,
                cwd=PROJECT_ROOT,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,  # Capture stderr for error reporting
                text=True,
                timeout=3600  # 1 hour timeout
            )
        else:
            # Let output stream directly to terminal (much faster!)
            result = subprocess.run(
                cmd,
                cwd=PROJECT_ROOT,
                # stdout and stderr default to None, which means inherit from parent
                timeout=3600  # 1 hour timeout
            )

        if not quiet:
            print(f"  " + "-" * 76)

        if result.returncode == 0:
            print(f"  ✓ Static solve completed successfully")
            return True
        else:
            print(f"  ✗ Static solve failed with return code {result.returncode}")
            if quiet and hasattr(result, 'stderr') and result.stderr:
                print(f"  Error output:\n{result.stderr[-1000:]}")  # Last 1000 chars
            return False

    except subprocess.TimeoutExpired:
        print(f"  ✗ Static solve timed out after 1 hour")
        return False
    except Exception as e:
        print(f"  ✗ Error running static solve: {e}")
        return False


def extract_values_from_csv(case_folder):
    """
    Extract energy values from the CSV output file.

    Args:
        case_folder (Path): Path to the case folder

    Returns:
        dict: Dictionary with elasticity.i variable names as keys and values as floats
              Returns None if CSV not found or extraction fails
    """
    csv_file = case_folder / CSV_FILENAME

    if not csv_file.exists():
        print(f"  ✗ ERROR: CSV file not found: {csv_file}")
        return None

    try:
        with open(csv_file, 'r') as f:
            reader = csv.DictReader(f)

            # Read all rows (we want the last row for steady-state values)
            rows = list(reader)

            if not rows:
                print(f"  ✗ ERROR: CSV file is empty")
                return None

            # Get the last row (final steady-state values)
            last_row = rows[-1]

            # Extract the required values using the column name mapping
            extracted_values = {}
            for elasticity_var_name, csv_column_name in CSV_COLUMN_NAMES.items():
                if csv_column_name in last_row:
                    extracted_values[elasticity_var_name] = float(last_row[csv_column_name])
                else:
                    print(f"  ✗ WARNING: Column '{csv_column_name}' not found in CSV")
                    print(f"  Available columns: {', '.join(last_row.keys())}")
                    return None

            return extracted_values

    except Exception as e:
        print(f"  ✗ ERROR reading CSV file: {e}")
        return None


def update_elasticity_file(case_folder, values):
    """
    Update the first three lines of elasticity.i with extracted values.

    Args:
        case_folder (Path): Path to the case folder
        values (dict): Dictionary with variable names and values

    Returns:
        bool: True if successful, False otherwise
    """
    elasticity_file = case_folder / "elasticity.i"

    if not elasticity_file.exists():
        print(f"  ✗ ERROR: elasticity.i not found in {case_folder.name}")
        return False

    try:
        # Read the existing file
        with open(elasticity_file, 'r') as f:
            lines = f.readlines()

        # Prepare the new first three lines (in the correct order)
        new_lines = []
        for var_name in ELASTICITY_VARIABLES:
            new_lines.append(f"{var_name} = {values[var_name]:.6e}\n")

        # Replace the first three lines
        lines[:3] = new_lines

        # Write back to file
        with open(elasticity_file, 'w') as f:
            f.writelines(lines)

        print(f"  ✓ Updated elasticity.i with new values:")
        for line in new_lines:
            print(f"    {line.strip()}")

        return True

    except Exception as e:
        print(f"  ✗ ERROR updating elasticity.i: {e}")
        return False


def process_case(case_folder, dry_run=False, skip_simulation=False, mpi_processes=8, quiet=False):
    """
    Process a single case: run static solve and update elasticity.i

    Args:
        case_folder (Path): Path to the case folder
        dry_run (bool): If True, only extract and display values without updating
        skip_simulation (bool): If True, skip running the simulation (assumes CSV exists)
        mpi_processes (int): Number of MPI processes to use
        quiet (bool): If True, suppress MOOSE output for speed

    Returns:
        bool: True if successful, False otherwise
    """
    print("\n" + "=" * 80)
    print(f"Processing: {case_folder.name}")
    print("=" * 80)

    # Step 1: Run static solve (unless skipped)
    if not skip_simulation:
        success = run_static_solve(case_folder, mpi_processes, quiet)
        if not success:
            print(f"  ✗ Failed to run static solve for {case_folder.name}")
            return False
    else:
        print(f"  ⊙ Skipping simulation (using existing CSV)")

    # Step 2: Extract values from CSV
    print(f"\n  Extracting values from CSV...")
    values = extract_values_from_csv(case_folder)

    if values is None:
        print(f"  ✗ Failed to extract values from CSV")
        return False

    print(f"  ✓ Extracted values:")
    for var_name, value in values.items():
        print(f"    {var_name} = {value:.6e}")

    # Step 3: Update elasticity.i (unless dry run)
    if dry_run:
        print(f"\n  ⊙ DRY RUN: Would update elasticity.i with above values")
        return True
    else:
        print(f"\n  Updating elasticity.i...")
        success = update_elasticity_file(case_folder, values)
        return success


# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Run static solve and update elasticity files for parametric cases'
    )

    # Mutually exclusive group for case selection
    case_group = parser.add_mutually_exclusive_group(required=True)
    case_group.add_argument('--case', type=str,
                           help='Process a specific case folder (e.g., case_cf1_domain1x)')
    case_group.add_argument('--pattern', type=str,
                           help='Process all cases matching pattern (e.g., "case_cf*")')
    case_group.add_argument('--all', action='store_true',
                           help='Process all case folders')

    # Optional arguments
    parser.add_argument('--dry-run', action='store_true',
                       help='Extract values but do not update elasticity.i')
    parser.add_argument('--skip-simulation', action='store_true',
                       help='Skip running simulation (use existing CSV files)')
    parser.add_argument('--quiet', '-q', action='store_true',
                       help='Suppress MOOSE output for faster execution')
    parser.add_argument('--np', type=int, default=8,
                       help=f'Number of MPI processes (default: 8)')

    args = parser.parse_args()

    # Use the specified number of MPI processes
    mpi_processes = args.np

    # Print header
    print("=" * 80)
    print("Static Solve Runner and Elasticity File Updater")
    print("=" * 80)
    print(f"\nScript directory: {SCRIPT_DIR}")
    print(f"Project root: {PROJECT_ROOT}")
    print(f"Executable: {EXECUTABLE}")
    print(f"MPI processes: {mpi_processes}")

    # Check if executable exists
    executable_path = PROJECT_ROOT / EXECUTABLE.lstrip('./')
    if not executable_path.exists():
        print(f"\n✗ ERROR: Executable not found: {executable_path}")
        print(f"  Please ensure farms-opt is built and located at the project root")
        sys.exit(1)

    # Determine which cases to process
    if args.case:
        case_folder = SCRIPT_DIR / args.case
        if not case_folder.is_dir():
            print(f"\n✗ ERROR: Case folder not found: {case_folder}")
            sys.exit(1)
        case_folders = [case_folder]
    elif args.pattern:
        case_folders = find_case_folders(args.pattern)
        if not case_folders:
            print(f"\n✗ ERROR: No case folders found matching pattern: {args.pattern}")
            sys.exit(1)
    else:  # --all
        case_folders = find_case_folders()
        if not case_folders:
            print(f"\n✗ ERROR: No case folders found")
            sys.exit(1)

    print(f"\nFound {len(case_folders)} case(s) to process:")
    for cf in case_folders:
        print(f"  - {cf.name}")

    if args.dry_run:
        print("\n⊙ DRY RUN MODE: Will not modify elasticity.i files")

    if args.skip_simulation:
        print("\n⊙ SKIP SIMULATION MODE: Will use existing CSV files")

    if args.quiet:
        print("\n⊙ QUIET MODE: MOOSE output will be suppressed for faster execution")

    # Process each case
    print("\n" + "=" * 80)
    print("Starting processing...")
    print("=" * 80)

    results = {}
    for case_folder in case_folders:
        success = process_case(case_folder,
                             dry_run=args.dry_run,
                             skip_simulation=args.skip_simulation,
                             mpi_processes=mpi_processes,
                             quiet=args.quiet)
        results[case_folder.name] = success

    # Print summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    successful = [name for name, success in results.items() if success]
    failed = [name for name, success in results.items() if not success]

    print(f"\nTotal cases: {len(results)}")
    print(f"Successful: {len(successful)}")
    print(f"Failed: {len(failed)}")

    if successful:
        print("\n✓ Successful cases:")
        for name in successful:
            print(f"  - {name}")

    if failed:
        print("\n✗ Failed cases:")
        for name in failed:
            print(f"  - {name}")

    print("\n" + "=" * 80)

    # Exit with appropriate code
    sys.exit(0 if not failed else 1)


if __name__ == "__main__":
    main()
