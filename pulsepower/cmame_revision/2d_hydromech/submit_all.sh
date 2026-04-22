#!/bin/bash
# =============================================================================
# Submit every Frontera job under 2d_hydromech/ in one shot.
#
# Runs on the Frontera login node after the inputs have been synced into
# scratch. Each hydromech case has two sbatch files:
#   submit_static.sbatch       (produces static_solve_out.e)
#   submit_elasticity.sbatch   (reads static_solve_out.e as initial condition)
# so this script submits them as a dependency chain (elasticity runs
# afterok:<static_jid>) unless --skip-static is given.
#
# Discovers every submit_static.sbatch at depth <= 4 and submits both scripts
# from the case directory so SLURM .o*/.e* logs and MOOSE output files land
# next to the input files.
#
# Cases discovered today:
#   benchmark/
#   permeability_formula/
#
# Usage:
#   ./submit_all.sh                 # submit everything (with confirmation)
#   ./submit_all.sh --yes           # submit without confirmation
#   ./submit_all.sh --dry-run       # show what would be submitted
#   ./submit_all.sh --skip-static   # only submit submit_elasticity.sbatch
#                                   # (use when static_solve_out.e already exists)
#   ./submit_all.sh benchmark       # only cases matching 'benchmark'
#   ./submit_all.sh permeability    # only cases matching 'permeability'
#
# Artifacts:
#   submitted_jobs.txt  — JOBID:case_path[:stage] for each successful
#                         submission. Used by monitor/cancel hints below.
# =============================================================================

set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

DRY_RUN=0
SKIP_CONFIRM=0
SKIP_STATIC=0
PATTERN=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)       DRY_RUN=1 ;;
        --yes|-y)        SKIP_CONFIRM=1 ;;
        --skip-static)   SKIP_STATIC=1 ;;
        -h|--help)
            sed -n '2,35p' "$0"
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

# Discover case directories by finding every submit_static.sbatch (each case
# has exactly one). -maxdepth 4 handles ./<case>/, ./<group>/<case>/, ...
mapfile -t STATIC_SCRIPTS < <(find "$SCRIPT_DIR" -maxdepth 4 -name 'submit_static.sbatch' -type f | sort)

if [[ -n "$PATTERN" ]]; then
    FILTERED=()
    for s in "${STATIC_SCRIPTS[@]}"; do
        if [[ "$s" == *"$PATTERN"* ]]; then
            FILTERED+=("$s")
        fi
    done
    STATIC_SCRIPTS=("${FILTERED[@]}")
fi

if [[ ${#STATIC_SCRIPTS[@]} -eq 0 ]]; then
    echo "ERROR: no submit_static.sbatch scripts matched ${PATTERN:-<everything>} under $SCRIPT_DIR" >&2
    exit 1
fi

# Each case must also have an elasticity sbatch; refuse to submit if it's missing.
for ss in "${STATIC_SCRIPTS[@]}"; do
    es="$(dirname "$ss")/submit_elasticity.sbatch"
    if [[ ! -f "$es" ]]; then
        echo "ERROR: missing $es (every case must have both sbatch files)" >&2
        exit 1
    fi
done

echo "Found ${#STATIC_SCRIPTS[@]} case(s):"
for s in "${STATIC_SCRIPTS[@]}"; do
    echo "  $(dirname "${s#$SCRIPT_DIR/}")"
done
if [[ $SKIP_STATIC -eq 1 ]]; then
    echo "Mode: --skip-static (submit elasticity only; no dependency)"
else
    echo "Mode: chain (static -> elasticity afterok:<static_jid>)"
fi

if [[ $DRY_RUN -eq 1 ]]; then
    echo "[dry-run] exiting without submitting."
    exit 0
fi

if [[ $SKIP_CONFIRM -eq 0 ]]; then
    echo
    n_jobs=$(( ${#STATIC_SCRIPTS[@]} * ( SKIP_STATIC == 1 ? 1 : 2 ) ))
    read -rp "Submit ${n_jobs} job(s)? [y/N] " REPLY
    case "$REPLY" in
        y|Y|yes|YES) ;;
        *) echo "cancelled."; exit 0 ;;
    esac
fi

SUBMITTED_IDS=()
FAILED=()

echo
echo "--- Submitting ---"
for static_script in "${STATIC_SCRIPTS[@]}"; do
    case_dir="$(dirname "$static_script")"
    case_name="${case_dir#$SCRIPT_DIR/}"
    pushd "$case_dir" > /dev/null

    static_jid=""
    if [[ $SKIP_STATIC -eq 0 ]]; then
        OUT="$(sbatch submit_static.sbatch 2>&1)" || true
        if [[ "$OUT" =~ Submitted\ batch\ job\ ([0-9]+) ]]; then
            static_jid="${BASH_REMATCH[1]}"
            printf '  [OK]   %-50s static      jobid=%s\n' "$case_name" "$static_jid"
            SUBMITTED_IDS+=("$static_jid:$case_name:static")
        else
            printf '  [FAIL] %-50s static      sbatch: %s\n' "$case_name" "$OUT"
            FAILED+=("$case_name:static")
            popd > /dev/null
            continue   # don't submit elasticity without a valid dependency
        fi
    fi

    if [[ -n "$static_jid" ]]; then
        OUT="$(sbatch --dependency=afterok:"$static_jid" submit_elasticity.sbatch 2>&1)" || true
    else
        OUT="$(sbatch submit_elasticity.sbatch 2>&1)" || true
    fi
    if [[ "$OUT" =~ Submitted\ batch\ job\ ([0-9]+) ]]; then
        JID="${BASH_REMATCH[1]}"
        dep_note="${static_jid:+(afterok:$static_jid)}"
        printf '  [OK]   %-50s elasticity  jobid=%s  %s\n' "$case_name" "$JID" "$dep_note"
        SUBMITTED_IDS+=("$JID:$case_name:elasticity")
    else
        printf '  [FAIL] %-50s elasticity  sbatch: %s\n' "$case_name" "$OUT"
        FAILED+=("$case_name:elasticity")
    fi

    popd > /dev/null
done

# Persist job IDs so monitor/cancel/sync scripts can use them.
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
