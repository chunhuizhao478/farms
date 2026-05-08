#!/bin/bash
# =============================================================================
# Submit every Frontera job under viscosity/ in one shot.
#
# Discovers every submit_elasticity.sbatch at depth <= 2 and submits each from
# its own case directory so the SLURM .o*/.e* logs and MOOSE output files land
# next to the input files.
#
# Cases:
#   pure_solid/   (no fluid, EM=0.100, 10 pulses) — drained-limit reference
#   mu_1em4/      (mu = 1e-4 Pa.s, EM=0.100, 10 pulses)
#   mu_1em3/      (mu = 1e-3 Pa.s, EM=0.100, 10 pulses) — water baseline
#   mu_1em1/      (mu = 1e-1 Pa.s, EM=0.100, 10 pulses) — undrained limit
#
# Each hydromech case reads ../../static_solve_out.e via SolutionUserObject;
# make sure that file exists on the cluster before submitting.
#
# Usage:
#   ./submit_all.sh                 # submit everything (with confirmation)
#   ./submit_all.sh --yes           # submit without confirmation
#   ./submit_all.sh --dry-run       # show what would be submitted
#   ./submit_all.sh mu_1em3         # only submit cases matching 'mu_1em3'
#
# Artifacts:
#   submitted_jobs.txt  — JOBID:case_path for each successful submission.
# =============================================================================

set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

DRY_RUN=0
SKIP_CONFIRM=0
PATTERN=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)       DRY_RUN=1 ;;
        --yes|-y)        SKIP_CONFIRM=1 ;;
        -h|--help)
            sed -n '2,24p' "$0"
            exit 0
            ;;
        --*)
            echo "Unknown option: $1" >&2
            exit 2
            ;;
        *)
            PATTERN="$1"
            ;;
    esac
    shift
done

# Discover every submit_elasticity.sbatch under each case directory.
SCRIPTS=()
while IFS= read -r line; do
    SCRIPTS+=("$line")
done < <(find "$SCRIPT_DIR" -maxdepth 2 \
              -name 'submit_elasticity.sbatch' -type f -print | sort)

if [[ -n "$PATTERN" ]]; then
    FILTERED=()
    for s in "${SCRIPTS[@]}"; do
        if [[ "$s" == *"$PATTERN"* ]]; then
            FILTERED+=("$s")
        fi
    done
    SCRIPTS=("${FILTERED[@]}")
fi

if [[ ${#SCRIPTS[@]} -eq 0 ]]; then
    echo "ERROR: no submit_elasticity.sbatch scripts matched ${PATTERN:-<everything>} under $SCRIPT_DIR" >&2
    exit 1
fi

echo "Found ${#SCRIPTS[@]} job script(s):"
for s in "${SCRIPTS[@]}"; do
    echo "  ${s#$SCRIPT_DIR/}"
done

if [[ $DRY_RUN -eq 1 ]]; then
    echo "[dry-run] exiting without submitting."
    exit 0
fi

if [[ $SKIP_CONFIRM -eq 0 ]]; then
    echo
    read -rp "Submit ${#SCRIPTS[@]} job(s)? [y/N] " REPLY
    case "$REPLY" in
        y|Y|yes|YES) ;;
        *) echo "cancelled."; exit 0 ;;
    esac
fi

SUBMITTED_IDS=()
FAILED=()

echo
echo "--- Submitting ---"
for script in "${SCRIPTS[@]}"; do
    case_dir="$(dirname "$script")"
    case_name="${case_dir#$SCRIPT_DIR/}"
    pushd "$case_dir" > /dev/null
    OUT="$(sbatch submit_elasticity.sbatch 2>&1)" || true
    popd > /dev/null
    if [[ "$OUT" =~ Submitted\ batch\ job\ ([0-9]+) ]]; then
        JID="${BASH_REMATCH[1]}"
        printf '  [OK]   %-40s  jobid=%s\n' "$case_name" "$JID"
        SUBMITTED_IDS+=("$JID:$case_name")
    else
        printf '  [FAIL] %-40s  sbatch: %s\n' "$case_name" "$OUT"
        FAILED+=("$case_name")
    fi
done

if [[ ${#SUBMITTED_IDS[@]} -gt 0 ]]; then
    printf '%s\n' "${SUBMITTED_IDS[@]}" > "$SCRIPT_DIR/submitted_jobs.txt"
fi

echo
echo "=== Summary ==="
echo "Submitted: ${#SUBMITTED_IDS[@]}"
echo "Failed:    ${#FAILED[@]}"
if [[ ${#SUBMITTED_IDS[@]} -gt 0 ]]; then
    echo
    echo "Job IDs recorded in submitted_jobs.txt"
    echo "Monitor:          squeue -u \$USER"
    echo "Cancel all above: scancel \$(awk -F: '{print \$1}' submitted_jobs.txt)"
fi
if [[ ${#FAILED[@]} -gt 0 ]]; then
    exit 1
fi
exit 0
