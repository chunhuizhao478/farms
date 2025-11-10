#!/usr/bin/env python3
"""
SLURM Job Submission Script for Parametric Studies

This script submits all generated SLURM job files for parametric study cases.
It can be used for both pure solid and hydro-mechanical (porousflow) simulations.

Usage:
    # Submit all jobs (with confirmation)
    python3 submit_all_jobs.py

    # Dry run (see what would be submitted without actually submitting)
    python3 submit_all_jobs.py --dry-run

    # Validate all required files exist before submission
    python3 submit_all_jobs.py --validate

    # Submit jobs matching a specific pattern
    python3 submit_all_jobs.py --pattern "case_cf*"

    # Submit a specific case
    python3 submit_all_jobs.py --case case_cf1_domain1x

    # Submit without confirmation prompt
    python3 submit_all_jobs.py --yes

    # Add delay between submissions
    python3 submit_all_jobs.py --delay 2

Author: Auto-generated
Date: 2025-11-09
"""

import argparse
import re
import subprocess
import sys
import time
from pathlib import Path

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Base directory (where this script is located)
BASE_DIR = Path(__file__).parent

# Default submit script name pattern
SUBMIT_SCRIPT_PATTERN = "submit_*.sh"

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


def find_submit_scripts(case_folder):
    """
    Find submit scripts in a case folder.

    Args:
        case_folder (Path): Path to the case folder

    Returns:
        list: List of Path objects for submit scripts
    """
    submit_scripts = sorted(case_folder.glob(SUBMIT_SCRIPT_PATTERN))
    return submit_scripts


def submit_job(script_path, dry_run=False):
    """
    Submit a SLURM job using sbatch.

    Args:
        script_path (Path): Path to the submit script
        dry_run (bool): If True, only print what would be done

    Returns:
        tuple: (success, job_id or error_message)
    """
    if dry_run:
        return True, "DRY-RUN"

    try:
        # Change to the directory containing the script
        original_dir = Path.cwd()
        script_dir = script_path.parent

        # Submit the job from its directory
        result = subprocess.run(
            ['sbatch', script_path.name],
            cwd=script_dir,
            capture_output=True,
            text=True,
            check=True
        )

        # Parse job ID from output (typically: "Submitted batch job 12345")
        output = result.stdout.strip()
        if "Submitted batch job" in output:
            job_id = output.split()[-1]
            return True, job_id
        else:
            return True, output

    except subprocess.CalledProcessError as e:
        return False, e.stderr.strip()
    except FileNotFoundError:
        return False, "sbatch command not found (not on HPC cluster?)"
    except Exception as e:
        return False, str(e)


def confirm_submission(num_jobs):
    """
    Ask user to confirm job submission.

    Args:
        num_jobs (int): Number of jobs to submit

    Returns:
        bool: True if user confirms, False otherwise
    """
    response = input(f"\nSubmit {num_jobs} job(s)? [y/N]: ").strip().lower()
    return response in ['y', 'yes']


def parse_submit_script(script_path):
    """
    Parse submit script to extract input file and executable paths.

    Args:
        script_path (Path): Path to the submit script

    Returns:
        dict: Dictionary with 'input_file' and 'executable' keys
    """
    info = {
        'input_file': None,
        'executable': None,
        'working_dir': None
    }

    try:
        with open(script_path, 'r') as f:
            content = f.read()

        # Find the ibrun/mpirun command with input file
        # Pattern: ibrun <executable> -i <input_file>
        run_match = re.search(r'(?:ibrun|mpirun|mpiexec)\s+(\S+)\s+-i\s+(\S+)', content)
        if run_match:
            info['executable'] = run_match.group(1)
            info['input_file'] = run_match.group(2)

    except Exception as e:
        pass

    return info


def find_referenced_files(input_file_path):
    """
    Find files referenced in MOOSE input file.

    Args:
        input_file_path (Path): Path to the input file

    Returns:
        list: List of referenced file paths
    """
    referenced_files = []

    if not input_file_path.exists():
        return referenced_files

    try:
        with open(input_file_path, 'r') as f:
            content = f.read()

        # Common patterns in MOOSE input files
        patterns = [
            r'file\s*=\s*["\']?([^"\'\s]+)',           # file = 'path'
            r'mesh_file\s*=\s*["\']?([^"\'\s]+)',      # mesh_file = 'path'
            r'data_file\s*=\s*["\']?([^"\'\s]+)',      # data_file = 'path'
            r'exodus\s*=\s*["\']?([^"\'\s]+)',         # exodus = 'path'
            r'file_base\s*=\s*["\']?([^"\'\s]+)',      # file_base = 'path'
        ]

        for pattern in patterns:
            matches = re.findall(pattern, content, re.IGNORECASE)
            for match in matches:
                # Remove quotes if present
                match = match.strip('"').strip("'")
                if match and not match.startswith('${'):  # Skip variables
                    referenced_files.append(match)

    except Exception as e:
        pass

    return list(set(referenced_files))  # Remove duplicates


def validate_job_files(case_folder, script_path, verbose=True):
    """
    Validate that all required files exist for a job.

    Args:
        case_folder (Path): Path to the case folder
        script_path (Path): Path to the submit script
        verbose (bool): If True, print detailed validation info

    Returns:
        tuple: (is_valid, missing_files, warnings)
    """
    missing_files = []
    warnings = []

    # Parse submit script
    script_info = parse_submit_script(script_path)

    # Check input file
    if script_info['input_file']:
        input_file_path = Path(script_info['input_file'])

        # For absolute paths (HPC paths), check the local case directory instead
        if input_file_path.is_absolute():
            # Extract just the filename and check in local case directory
            local_input_path = case_folder / input_file_path.name
        else:
            # Relative path - check relative to case folder
            local_input_path = case_folder / input_file_path.name

        if not local_input_path.exists():
            missing_files.append(('Input file', input_file_path.name))
        else:
            # Check files referenced in input file (use local path)
            referenced = find_referenced_files(local_input_path)
            for ref_file in referenced:
                ref_path = Path(ref_file)

                # Try multiple possible locations
                possible_paths = [
                    ref_path,  # Absolute path
                    case_folder / ref_path,  # Relative to case folder
                    case_folder / ref_path.name,  # Just filename in case folder
                    case_folder.parent / ref_path,  # Relative to parent
                    case_folder.parent.parent / ref_path,  # Two levels up
                ]

                file_found = False
                for possible_path in possible_paths:
                    if possible_path.exists():
                        file_found = True
                        break

                if not file_found:
                    # Only warn about mesh files, as others might be generated
                    if any(ext in ref_file.lower() for ext in ['.msh', '.e', '.exo', '.mesh']):
                        missing_files.append(('Referenced mesh file', ref_file))
                    else:
                        warnings.append(f"Referenced file not found (may be generated): {ref_file}")
    else:
        warnings.append("Could not parse input file path from submit script")

    # Check executable (might be in PATH or relative)
    if script_info['executable']:
        exe_name = script_info['executable']
        exe_path = case_folder / exe_name

        # Don't check if it's just a command name (will be in PATH)
        if '/' in exe_name and not exe_path.exists():
            # Check parent directories
            parent_exe = case_folder.parent / exe_name
            if not parent_exe.exists():
                warnings.append(f"Executable not found locally (may be in PATH): {exe_name}")

    is_valid = len(missing_files) == 0

    if verbose:
        if is_valid:
            status = "✓"
        else:
            status = "✗"

        print(f"  {status} [{case_folder.name}] {script_path.name}")

        if missing_files:
            for file_type, file_path in missing_files:
                print(f"      ✗ Missing {file_type}: {file_path}")

        if warnings and verbose:
            for warning in warnings:
                print(f"      ⚠ {warning}")

    return is_valid, missing_files, warnings


# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Submit SLURM jobs for all parametric cases',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Submit all jobs with confirmation
  python3 submit_all_jobs.py

  # Dry run to see what would be submitted
  python3 submit_all_jobs.py --dry-run

  # Submit only cases matching pattern
  python3 submit_all_jobs.py --pattern "case_cf*"

  # Submit without confirmation
  python3 submit_all_jobs.py --yes

  # Add 2 second delay between submissions
  python3 submit_all_jobs.py --delay 2
        """
    )

    # Case selection
    case_group = parser.add_mutually_exclusive_group()
    case_group.add_argument('--case', type=str,
                           help='Submit jobs for a specific case folder')
    case_group.add_argument('--pattern', type=str,
                           help='Submit jobs for all cases matching pattern (e.g., "case_cf*")')

    # Submission options
    parser.add_argument('--dry-run', action='store_true',
                       help='Show what would be submitted without actually submitting')
    parser.add_argument('--validate', '--check', action='store_true',
                       help='Validate that all required files exist before submission')
    parser.add_argument('--yes', '-y', action='store_true',
                       help='Skip confirmation prompt')
    parser.add_argument('--delay', type=float, default=0,
                       help='Delay in seconds between job submissions (default: 0)')
    parser.add_argument('--script-name', type=str,
                       help='Submit only scripts matching this name (e.g., submit_elasticity.sh)')

    args = parser.parse_args()

    # Print header
    print("=" * 80)
    print("SLURM Job Submission Script")
    print("=" * 80)
    print(f"\nBase directory: {BASE_DIR}")
    if args.dry_run:
        print("\n*** DRY RUN MODE - No jobs will actually be submitted ***")

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

    # Find all submit scripts
    jobs_to_submit = []
    for case_folder in case_folders:
        submit_scripts = find_submit_scripts(case_folder)

        # Filter by script name if specified
        if args.script_name:
            submit_scripts = [s for s in submit_scripts if s.name == args.script_name]

        for script in submit_scripts:
            jobs_to_submit.append((case_folder, script))

    if not jobs_to_submit:
        print(f"\n✗ ERROR: No submit scripts found")
        if args.script_name:
            print(f"   (looking for: {args.script_name})")
        sys.exit(1)

    # Display jobs to be submitted
    print(f"\nFound {len(jobs_to_submit)} job(s) to submit:")
    print()
    for case_folder, script in jobs_to_submit:
        print(f"  [{case_folder.name}] {script.name}")

    # Validate files if requested
    if args.validate:
        print("\n" + "=" * 80)
        print("Validating required files...")
        print("=" * 80)
        print()

        all_valid = True
        validation_results = []

        for case_folder, script in jobs_to_submit:
            is_valid, missing, warnings = validate_job_files(case_folder, script, verbose=True)
            validation_results.append((case_folder, script, is_valid, missing, warnings))
            if not is_valid:
                all_valid = False

        # Print validation summary
        print("\n" + "=" * 80)
        print("Validation Summary")
        print("=" * 80)

        valid_count = sum(1 for _, _, is_valid, _, _ in validation_results if is_valid)
        invalid_count = len(validation_results) - valid_count

        print(f"\nTotal jobs: {len(validation_results)}")
        print(f"Valid: {valid_count}")
        print(f"Invalid (missing files): {invalid_count}")

        if not all_valid:
            print("\n✗ Validation FAILED - Some required files are missing")
            print("\nJobs with missing files:")
            for case_folder, script, is_valid, missing, _ in validation_results:
                if not is_valid:
                    print(f"  [{case_folder.name}] {script.name}")
                    for file_type, file_path in missing:
                        print(f"    - {file_type}: {file_path}")
            sys.exit(1)
        else:
            print("\n✓ Validation PASSED - All required files are present")

        # If only validating, exit here
        if args.dry_run or (not args.yes and args.validate):
            print("\nValidation complete. Use --yes to proceed with submission.")
            sys.exit(0)

    # Confirm submission unless --yes or --dry-run
    if not args.yes and not args.dry_run:
        if not confirm_submission(len(jobs_to_submit)):
            print("\nSubmission cancelled.")
            sys.exit(0)

    # Submit jobs
    print("\n" + "=" * 80)
    print("Submitting jobs...")
    print("=" * 80)
    print()

    submitted = []
    failed = []

    for i, (case_folder, script) in enumerate(jobs_to_submit, 1):
        # Add delay between submissions (except for first job)
        if args.delay > 0 and i > 1:
            time.sleep(args.delay)

        # Submit the job
        success, result = submit_job(script, dry_run=args.dry_run)

        if success:
            submitted.append((case_folder, script, result))
            status_icon = "✓" if not args.dry_run else "○"
            print(f"  {status_icon} [{case_folder.name}] {script.name} → Job ID: {result}")
        else:
            failed.append((case_folder, script, result))
            print(f"  ✗ [{case_folder.name}] {script.name} → ERROR: {result}")

    # Print summary
    print("\n" + "=" * 80)
    print("Submission Summary")
    print("=" * 80)
    print(f"\nTotal jobs: {len(jobs_to_submit)}")
    print(f"Successfully submitted: {len(submitted)}")
    print(f"Failed: {len(failed)}")

    if failed:
        print("\nFailed submissions:")
        for case_folder, script, error in failed:
            print(f"  ✗ [{case_folder.name}] {script.name}")
            print(f"    Error: {error}")

    if submitted and not args.dry_run:
        print("\nMonitoring commands:")
        print("  # Check job status")
        print("  squeue -u $USER")
        print()
        print("  # Check detailed job info")
        print(f"  scontrol show job {submitted[0][2]}")
        print()
        print("  # Cancel all jobs")
        print("  scancel -u $USER")
        print()
        print("  # Cancel specific job")
        print(f"  scancel {submitted[0][2]}")

    # Exit with error code if any submissions failed
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
