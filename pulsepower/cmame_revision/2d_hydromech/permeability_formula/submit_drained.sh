#!/bin/bash
# =============================================================================
# Submit the three drained-condition jobs (3 pulses, end_time=3e-5) in one shot.
#
# Run on the Frontera login node after the inputs are synced to scratch.
# Each case ships submit_elasticity_drained.sbatch, which reads the shared
# ../static_solve_out.e (permeability_formula/static_solve_out.e) as its IC.
#
# Cases:
#   undrained/           standard crack normal (damage_gradient),     DRAINED BC
#   strained-based/      principal_strain normal (Liu 2024),          DRAINED BC
#   regularized_normal/  regularized damage_gradient (eps = 1e-8),    DRAINED BC
#
# Usage:
#   ./submit_drained.sh            # submit (with confirmation)
#   ./submit_drained.sh --yes      # submit without confirmation
#   ./submit_drained.sh --dry-run  # show what would be submitted
#
# Records JOBID:case for each submission to submitted_jobs_drained.txt.
# =============================================================================

set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASES=(undrained strained-based regularized_normal)
JOB_SCRIPT=submit_elasticity_drained.sbatch
STATIC_REF=static_solve_out.e   # shared, one level up from each case dir

DRY_RUN=0
SKIP_CONFIRM=0
while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run) DRY_RUN=1 ;;
        --yes|-y)  SKIP_CONFIRM=1 ;;
        -h|--help) sed -n '2,20p' "$0"; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
    shift
done

# Every case directory must ship its drained sbatch.
MISSING=0
for c in "${CASES[@]}"; do
    if [[ ! -f "$SCRIPT_DIR/$c/$JOB_SCRIPT" ]]; then
        echo "ERROR: missing $c/$JOB_SCRIPT" >&2
        MISSING=1
    fi
done
[[ $MISSING -eq 1 ]] && exit 1

# Shared static-solve IC sanity check (jobs read ../static_solve_out.e).
if [[ ! -f "$SCRIPT_DIR/$STATIC_REF" ]]; then
    echo "[WARN] $STATIC_REF not found in $SCRIPT_DIR; the jobs read ../$STATIC_REF and will fail until it is uploaded."
fi

echo "Drained jobs to submit (from $SCRIPT_DIR):"
for c in "${CASES[@]}"; do printf '  %s/%s\n' "$c" "$JOB_SCRIPT"; done

if [[ $DRY_RUN -eq 1 ]]; then
    echo "[dry-run] exiting without submitting."
    exit 0
fi

if [[ $SKIP_CONFIRM -eq 0 ]]; then
    echo
    read -rp "Submit ${#CASES[@]} job(s)? [y/N] " REPLY
    case "$REPLY" in
        y|Y|yes|YES) ;;
        *) echo "cancelled."; exit 0 ;;
    esac
fi

SUBMITTED=()
FAILED=()

echo
echo "--- Submitting ---"
for c in "${CASES[@]}"; do
    pushd "$SCRIPT_DIR/$c" > /dev/null
    OUT="$(sbatch "$JOB_SCRIPT" 2>&1)" || true
    if [[ "$OUT" =~ Submitted\ batch\ job\ ([0-9]+) ]]; then
        JID="${BASH_REMATCH[1]}"
        printf '  [OK]   %-20s jobid=%s\n' "$c" "$JID"
        SUBMITTED+=("$JID:$c:elasticity_drained")
    else
        printf '  [FAIL] %-20s sbatch: %s\n' "$c" "$OUT"
        FAILED+=("$c")
    fi
    popd > /dev/null
done

if [[ ${#SUBMITTED[@]} -gt 0 ]]; then
    printf '%s\n' "${SUBMITTED[@]}" > "$SCRIPT_DIR/submitted_jobs_drained.txt"
fi

echo
echo "=== Summary ==="
echo "Submitted: ${#SUBMITTED[@]}"
echo "Failed:    ${#FAILED[@]}"
if [[ ${#SUBMITTED[@]} -gt 0 ]]; then
    echo
    echo "Job IDs recorded in submitted_jobs_drained.txt"
    echo "Monitor:          squeue -u \$USER"
    echo "Cancel all above: scancel \$(awk -F: '{print \$1}' submitted_jobs_drained.txt)"
fi
[[ ${#FAILED[@]} -gt 0 ]] && exit 1
exit 0
