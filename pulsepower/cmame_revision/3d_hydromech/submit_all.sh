#!/bin/bash
# =============================================================================
# Submit every Frontera job under 3d_hydromech/ in one shot.
#
# Discovers every submit_elasticity.sbatch at depth <= 3 and submits each
# from its own case directory so the SLURM .o*/.e* logs and MOOSE output
# files land next to the input files.
#
# Cases discovered today:
#   em_25J/          (EM=0.025, undrained borehole)
#   em_50J/          (EM=0.050, undrained borehole)
#   em_5J_drained/   (EM=0.005, drained borehole: pp=func_tri_pulse on boundary 3)
#   em_10J_drained/  (EM=0.010, drained borehole: pp=func_tri_pulse on boundary 3)
#
# Usage:
#   ./submit_all.sh                 # submit everything (with confirmation)
#   ./submit_all.sh --yes           # submit without confirmation
#   ./submit_all.sh --dry-run       # show what would be submitted
#   ./submit_all.sh em_25J          # only submit cases matching 'em_25J'
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
            sed -n '2,22p' "$0"
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

# Discover every submit_elasticity.sbatch.  Skip the archive/ tree.
mapfile -t SCRIPTS < <(find "$SCRIPT_DIR" -maxdepth 3 -path "*/archive/*" -prune -o \
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
        printf '  [OK]   %-55s  jobid=%s\n' "$case_name" "$JID"
        SUBMITTED_IDS+=("$JID:$case_name")
    else
        printf '  [FAIL] %-55s  sbatch: %s\n' "$case_name" "$OUT"
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
