#!/bin/bash
# =============================================================================
# Submit every Frontera job in pulse_duration_magnitude/ in one shot.
#
# Discovers every submit_elasticity_*.sbatch in this directory and submits
# each from this directory so the SLURM .o*/.e* logs land here.
#
# Cases:
#   submit_elasticity_E1d25.sbatch              (EM=0.00125, dt_pulse=1e-5)
#   submit_elasticity_E1d25_pulse2em5.sbatch    (EM=0.00125, dt_pulse=2e-5)
#   submit_elasticity_E1d25_pulse4em5.sbatch    (EM=0.00125, dt_pulse=4e-5)
#   submit_elasticity_E2d50.sbatch              (EM=0.0025,  dt_pulse=1e-5)
#   submit_elasticity_E5d00.sbatch              (EM=0.005,   dt_pulse=1e-5)
#
# Each elasticity_*.i reads ../static_solve_out.e via SolutionUserObject;
# make sure that file exists on the cluster before submitting.
#
# Usage:
#   ./submit_all.sh                 # submit everything (with confirmation)
#   ./submit_all.sh --yes           # submit without confirmation
#   ./submit_all.sh --dry-run       # show what would be submitted
#   ./submit_all.sh E2d50           # only submit cases matching 'E2d50'
#
# Artifacts:
#   submitted_jobs.txt  — JOBID:case_name for each successful submission.
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
            sed -n '2,30p' "$0"
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

SCRIPTS=()
while IFS= read -r line; do
    SCRIPTS+=("$line")
done < <(find "$SCRIPT_DIR" -maxdepth 1 -name 'submit_elasticity_*.sbatch' -type f -print | sort)

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
    echo "ERROR: no submit_elasticity_*.sbatch matched ${PATTERN:-<everything>} under $SCRIPT_DIR" >&2
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
    case_name="$(basename "$script" .sbatch)"
    case_name="${case_name#submit_elasticity_}"
    pushd "$SCRIPT_DIR" > /dev/null
    OUT="$(sbatch "$(basename "$script")" 2>&1)" || true
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
