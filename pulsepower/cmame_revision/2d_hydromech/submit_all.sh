#!/bin/bash
# =============================================================================
# Submit every Frontera job under 2d_hydromech/ in one shot.
#
# Runs on the Frontera login node after the inputs have been synced into
# scratch. Each hydromech case has two sbatch files:
#   submit_static.sbatch       (produces static_solve_out.e)
#   submit_elasticity.sbatch   (reads static_solve_out.e as initial condition)
# This script submits ONLY submit_elasticity.sbatch. The static solve is
# expected to have been run elsewhere and its static_solve_out.e uploaded
# into the case directory before the elasticity job starts.
#
# Use --with-static to also submit the static job first and chain the
# elasticity job with --dependency=afterok:<static_jid>.
#
# Discovers every submit_elasticity.sbatch at depth <= 5 and submits it
# from the case directory so SLURM .o*/.e* logs and MOOSE output files land
# next to the input files.
#
# Cases discovered today (only under permeability_formula/):
#   permeability_formula/                                 (original single case)
#   permeability_formula/undrained/em_{0p005..0p100}/     (EM sweep, undrained)
#   permeability_formula/drained/em_{0p005..0p100}/       (EM sweep, drained)
#   permeability_formula/order_test/em_0p005/             (2nd-order disp + 1st-order pp)
# All sub-cases share ../../static_solve_out.e (must be uploaded to the
# parent permeability_formula/ on Frontera before submitting). Other
# subdirectories of 2d_hydromech (e.g. benchmark/) are NOT submitted by
# this script.
#
# Usage:
#   ./submit_all.sh                 # submit elasticity only (with confirmation)
#   ./submit_all.sh --yes           # submit without confirmation
#   ./submit_all.sh --dry-run       # show what would be submitted
#   ./submit_all.sh --with-static   # also submit static (chain via afterok)
#   ./submit_all.sh undrained       # only undrained EM sweep
#   ./submit_all.sh drained         # only drained EM sweep
#   ./submit_all.sh em_0p040        # only the EM=0.040 cases (both drained and undrained)
#   ./submit_all.sh order_test      # only the mixed-order test
#
# Artifacts:
#   submitted_jobs.txt  — JOBID:case_path:stage for each successful
#                         submission. Used by monitor/cancel hints below.
# =============================================================================

set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Restrict discovery to the permeability_formula tree only. Other subdirs
# of 2d_hydromech (e.g. benchmark/) are not submitted by this script.
SEARCH_ROOT="$SCRIPT_DIR/permeability_formula"
if [[ ! -d "$SEARCH_ROOT" ]]; then
    echo "ERROR: search root not found: $SEARCH_ROOT" >&2
    exit 1
fi

DRY_RUN=0
SKIP_CONFIRM=0
WITH_STATIC=0
PATTERN=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)       DRY_RUN=1 ;;
        --yes|-y)        SKIP_CONFIRM=1 ;;
        --with-static)   WITH_STATIC=1 ;;
        -h|--help)
            sed -n '2,37p' "$0"
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

# Discover case directories by finding every submit_elasticity.sbatch
# (each case has exactly one). -maxdepth 5 handles ./<case>/, ./<group>/<case>/, ...
# Portable array fill for bash 3.2+ (mapfile is bash 4+ only).
ELAST_SCRIPTS=()
while IFS= read -r s; do
    ELAST_SCRIPTS+=("$s")
done < <(find "$SEARCH_ROOT" -maxdepth 5 -name 'submit_elasticity.sbatch' -type f | sort)

if [[ -n "$PATTERN" ]]; then
    FILTERED=()
    for s in "${ELAST_SCRIPTS[@]}"; do
        if [[ "$s" == *"$PATTERN"* ]]; then
            FILTERED+=("$s")
        fi
    done
    ELAST_SCRIPTS=("${FILTERED[@]}")
fi

if [[ ${#ELAST_SCRIPTS[@]} -eq 0 ]]; then
    echo "ERROR: no submit_elasticity.sbatch scripts matched ${PATTERN:-<everything>} under $SCRIPT_DIR" >&2
    exit 1
fi

# If --with-static, every case must also have submit_static.sbatch.
if [[ $WITH_STATIC -eq 1 ]]; then
    for es in "${ELAST_SCRIPTS[@]}"; do
        ss="$(dirname "$es")/submit_static.sbatch"
        if [[ ! -f "$ss" ]]; then
            echo "ERROR: --with-static requested but missing $ss" >&2
            exit 1
        fi
    done
fi

echo "Found ${#ELAST_SCRIPTS[@]} case(s):"
for s in "${ELAST_SCRIPTS[@]}"; do
    echo "  $(dirname "${s#$SCRIPT_DIR/}")"
done
if [[ $WITH_STATIC -eq 1 ]]; then
    echo "Mode: chain (static -> elasticity afterok:<static_jid>)"
else
    echo "Mode: elasticity only (static_solve_out.e expected in each case dir)"
fi

if [[ $DRY_RUN -eq 1 ]]; then
    echo "[dry-run] exiting without submitting."
    exit 0
fi

if [[ $SKIP_CONFIRM -eq 0 ]]; then
    echo
    n_jobs=$(( ${#ELAST_SCRIPTS[@]} * ( WITH_STATIC == 1 ? 2 : 1 ) ))
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
for elast_script in "${ELAST_SCRIPTS[@]}"; do
    case_dir="$(dirname "$elast_script")"
    case_name="${case_dir#$SCRIPT_DIR/}"
    pushd "$case_dir" > /dev/null

    static_jid=""
    if [[ $WITH_STATIC -eq 1 ]]; then
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
    else
        # Sanity check: resolve the static-solve output path referenced by
        # this case's elasticity input (SolutionUserObject 'mesh = ...e').
        # Some cases use ./static_solve_out.e, others use ../../static_solve_out.e
        # (shared parent). Warn only if the referenced file isn't actually there.
        static_ref=""
        for input_i in elasticity*.i; do
            [[ -f "$input_i" ]] || continue
            static_ref=$(grep -Eo 'mesh = [^[:space:]]+\.e' "$input_i" | head -1 | sed 's|^mesh = ||')
            [[ -n "$static_ref" ]] && break
        done
        if [[ -n "$static_ref" ]]; then
            if [[ ! -f "$static_ref" ]]; then
                echo "  [WARN] $case_name: referenced static-solve output '$static_ref' not found (upload it before the job starts)."
            fi
        elif [[ ! -f "static_solve_out.e" ]]; then
            # No SolutionUserObject reference parsed; fall back to the
            # historical check against ./static_solve_out.e.
            echo "  [WARN] $case_name: static_solve_out.e not found; elasticity job will fail until you upload it."
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
