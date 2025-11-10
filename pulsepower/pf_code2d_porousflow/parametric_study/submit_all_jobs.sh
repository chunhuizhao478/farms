#!/bin/bash
################################################################################
# SLURM Job Submission Script (Bash Version)
#
# Simple bash script to submit all SLURM jobs for parametric studies.
# This script finds all submit_*.sh scripts in case_* directories and submits them.
#
# Usage:
#   ./submit_all_jobs.sh                    # Submit all jobs
#   ./submit_all_jobs.sh --validate         # Validate files without submitting
#   ./submit_all_jobs.sh case_cf*           # Submit jobs matching pattern
#   ./submit_all_jobs.sh case_cf1_domain1x  # Submit specific case only
#
# Options:
#   --validate, --check     Check if all required files exist
#   --dry-run              Show what would be submitted without submitting
#   --yes, -y              Skip confirmation prompt
#
# Author: Auto-generated
# Date: 2025-11-09
################################################################################

set -e  # Exit on error

# Colors for output (optional, comment out if not supported)
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Default options
VALIDATE_ONLY=0
DRY_RUN=0
SKIP_CONFIRM=0
PATTERN="case_*"

# Delay between submissions (in seconds)
DELAY=0

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --validate|--check)
            VALIDATE_ONLY=1
            shift
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        --yes|-y)
            SKIP_CONFIRM=1
            shift
            ;;
        --help|-h)
            echo "Usage: $0 [OPTIONS] [PATTERN]"
            echo ""
            echo "Options:"
            echo "  --validate, --check    Validate that all required files exist"
            echo "  --dry-run             Show what would be submitted without submitting"
            echo "  --yes, -y             Skip confirmation prompt"
            echo "  --help, -h            Show this help message"
            echo ""
            echo "Examples:"
            echo "  $0                       # Submit all jobs"
            echo "  $0 --validate            # Check all files exist"
            echo "  $0 case_cf*              # Submit jobs matching pattern"
            echo "  $0 --yes case_cf1_*      # Submit without confirmation"
            exit 0
            ;;
        *)
            PATTERN="$1"
            shift
            ;;
    esac
done

################################################################################
# Helper Functions
################################################################################

print_header() {
    echo "================================================================================"
    echo "$1"
    echo "================================================================================"
}

print_section() {
    echo ""
    echo "--------------------------------------------------------------------------------"
    echo "$1"
    echo "--------------------------------------------------------------------------------"
}

validate_job_files() {
    local case_dir="$1"
    local script_path="$2"
    local case_name=$(basename "$case_dir")
    local script_name=$(basename "$script_path")

    local is_valid=1
    local missing_files=()

    # Extract input file and executable from submit script
    local input_file=$(grep -oP '(?<=\s-i\s)\S+' "$script_path" 2>/dev/null | head -1)
    local executable=$(grep -oP '(?<=ibrun|mpirun|mpiexec\s)\S+(?=\s+-i)' "$script_path" 2>/dev/null | head -1)

    # Check if input file exists
    if [ -n "$input_file" ]; then
        # Try to find input file (could be absolute or relative path)
        local input_basename=$(basename "$input_file")
        local input_path=""

        if [ -f "$input_file" ]; then
            input_path="$input_file"
        elif [ -f "$case_dir/$input_basename" ]; then
            input_path="$case_dir/$input_basename"
        fi

        if [ -z "$input_path" ]; then
            is_valid=0
            missing_files+=("Input file: $input_file")
        else
            # Check for referenced mesh files in the input file
            local mesh_files=$(grep -oP '(?<=file\s=\s)["\047]?([^"\047\s]+\.(?:msh|e|exo|mesh))["\047]?' "$input_path" 2>/dev/null | tr -d '"' | tr -d "'")

            for mesh_file in $mesh_files; do
                local mesh_basename=$(basename "$mesh_file")
                local mesh_found=0

                # Check multiple possible locations
                for possible_dir in "$case_dir" "$case_dir/.." "$case_dir/../.." "$SCRIPT_DIR"; do
                    if [ -f "$possible_dir/$mesh_basename" ] || [ -f "$possible_dir/$mesh_file" ]; then
                        mesh_found=1
                        break
                    fi
                done

                if [ $mesh_found -eq 0 ]; then
                    is_valid=0
                    missing_files+=("Mesh file: $mesh_file")
                fi
            done
        fi
    fi

    # Print validation result
    if [ $is_valid -eq 1 ]; then
        echo -e "  ${GREEN}✓${NC} [$case_name] $script_name"
    else
        echo -e "  ${RED}✗${NC} [$case_name] $script_name"
        for missing in "${missing_files[@]}"; do
            echo -e "      ${RED}✗${NC} Missing $missing"
        done
    fi

    return $is_valid
}

################################################################################
# Main Script
################################################################################

print_header "SLURM Job Submission Script"

echo ""
echo "Base directory: $SCRIPT_DIR"
echo "Case pattern: $PATTERN"

if [ $VALIDATE_ONLY -eq 1 ]; then
    echo ""
    echo -e "${YELLOW}*** VALIDATION MODE - Checking files only ***${NC}"
elif [ $DRY_RUN -eq 1 ]; then
    echo ""
    echo -e "${YELLOW}*** DRY RUN MODE - No jobs will be submitted ***${NC}"
fi

# Change to the script directory
cd "$SCRIPT_DIR"

# Find all matching case directories
CASE_DIRS=($(ls -d $PATTERN 2>/dev/null | sort))

if [ ${#CASE_DIRS[@]} -eq 0 ]; then
    echo -e "${RED}✗ ERROR: No case directories found matching pattern: $PATTERN${NC}"
    exit 1
fi

echo ""
echo "Found ${#CASE_DIRS[@]} case(s) to process:"
for case_dir in "${CASE_DIRS[@]}"; do
    echo "  - $case_dir"
done

# Find all jobs to submit
JOBS_TO_SUBMIT=()
for case_dir in "${CASE_DIRS[@]}"; do
    if [ -d "$case_dir" ]; then
        # Find all submit scripts in this case directory
        for script in "$case_dir"/submit_*.sh; do
            if [ -f "$script" ]; then
                JOBS_TO_SUBMIT+=("$script")
            fi
        done
    fi
done

if [ ${#JOBS_TO_SUBMIT[@]} -eq 0 ]; then
    echo -e "${RED}✗ ERROR: No submit scripts found in case directories${NC}"
    exit 1
fi

print_section "Jobs to Submit"
for job in "${JOBS_TO_SUBMIT[@]}"; do
    case_name=$(dirname "$job")
    script_name=$(basename "$job")
    echo "  [$case_name] $script_name"
done

# Validate files if requested
if [ $VALIDATE_ONLY -eq 1 ]; then
    print_section "Validating Required Files"

    VALID_COUNT=0
    INVALID_COUNT=0

    for job in "${JOBS_TO_SUBMIT[@]}"; do
        case_dir=$(dirname "$job")

        if validate_job_files "$case_dir" "$job"; then
            ((VALID_COUNT++))
        else
            ((INVALID_COUNT++))
        fi
    done

    print_header "Validation Summary"
    echo ""
    echo "Total jobs: ${#JOBS_TO_SUBMIT[@]}"
    echo -e "Valid: ${GREEN}$VALID_COUNT${NC}"
    if [ $INVALID_COUNT -gt 0 ]; then
        echo -e "Invalid (missing files): ${RED}$INVALID_COUNT${NC}"
        echo ""
        echo -e "${RED}✗ Validation FAILED - Some required files are missing${NC}"
        exit 1
    else
        echo "Invalid: 0"
        echo ""
        echo -e "${GREEN}✓ Validation PASSED - All required files are present${NC}"
    fi

    print_header "Done"
    exit 0
fi

# Ask for confirmation unless --yes or --dry-run
if [ $SKIP_CONFIRM -eq 0 ] && [ $DRY_RUN -eq 0 ]; then
    echo ""
    read -p "Submit ${#JOBS_TO_SUBMIT[@]} job(s)? [y/N]: " -n 1 -r
    echo ""

    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Submission cancelled."
        exit 0
    fi
fi

print_section "Submitting Jobs"

SUBMITTED=0
FAILED=0
JOB_IDS=()

for job in "${JOBS_TO_SUBMIT[@]}"; do
    case_dir=$(dirname "$job")
    script_name=$(basename "$job")
    case_name=$(basename "$case_dir")

    if [ $DRY_RUN -eq 1 ]; then
        # Dry run - just show what would be submitted
        echo -e "  ${YELLOW}○${NC} [$case_name] $script_name → DRY-RUN"
        ((SUBMITTED++))
    else
        # Change to case directory and submit
        cd "$SCRIPT_DIR/$case_dir"

        # Submit the job
        OUTPUT=$(sbatch "$script_name" 2>&1)
        EXIT_CODE=$?

        if [ $EXIT_CODE -eq 0 ]; then
            # Extract job ID (format: "Submitted batch job 12345")
            JOB_ID=$(echo "$OUTPUT" | grep -oP "Submitted batch job \K\d+" || echo "$OUTPUT")
            echo -e "  ${GREEN}✓${NC} [$case_name] $script_name → Job ID: $JOB_ID"
            JOB_IDS+=("$JOB_ID")
            ((SUBMITTED++))
        else
            echo -e "  ${RED}✗${NC} [$case_name] $script_name → ERROR: $OUTPUT"
            ((FAILED++))
        fi

        # Add delay between submissions if specified
        if [ $DELAY -gt 0 ] && [ $SUBMITTED -lt ${#JOBS_TO_SUBMIT[@]} ]; then
            sleep $DELAY
        fi

        # Return to script directory
        cd "$SCRIPT_DIR"
    fi
done

# Print summary
print_header "Submission Summary"
echo ""
echo "Total jobs: ${#JOBS_TO_SUBMIT[@]}"
echo -e "Successfully submitted: ${GREEN}$SUBMITTED${NC}"
if [ $FAILED -gt 0 ]; then
    echo -e "Failed: ${RED}$FAILED${NC}"
else
    echo "Failed: 0"
fi

if [ $SUBMITTED -gt 0 ] && [ $DRY_RUN -eq 0 ]; then
    echo ""
    echo "Submitted Job IDs: ${JOB_IDS[@]}"

    print_section "Monitoring Commands"
    echo "  # Check job status"
    echo "  squeue -u \$USER"
    echo ""
    echo "  # Check detailed info for a specific job"
    echo "  scontrol show job ${JOB_IDS[0]}"
    echo ""
    echo "  # Watch job queue (updates every 2 seconds)"
    echo "  watch -n 2 'squeue -u \$USER'"
    echo ""
    echo "  # Check output files"
    echo "  tail -f case_*/*.o*"
    echo ""
    echo "  # Cancel all your jobs"
    echo "  scancel -u \$USER"
    echo ""
    echo "  # Cancel a specific job"
    echo "  scancel ${JOB_IDS[0]}"
    echo ""
fi

print_header "Done"

# Exit with error if any submissions failed
if [ $FAILED -gt 0 ]; then
    exit 1
fi

exit 0
