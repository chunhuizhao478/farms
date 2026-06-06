#!/bin/bash
# =============================================================================
# Submit every Frontera job under 2d_hydromech/ in one shot.
#
# Runs on the Frontera login node after the inputs have been synced into
# scratch. Each case ships elasticity job(s) that read a static-solve solution
# as their initial condition (SolutionUserObject). Cases in the 5-case
# comparison set ship TWO elasticity jobs (base + refined mesh):
#   submit_elasticity.sbatch          (reads ../static_solve_out.e,         base mesh)
#   submit_elasticity_refined.sbatch  (reads ../static_solve_refined_out.e, refined mesh)
# This script submits BOTH of those (older sweep cases ship only the base one).
# The static solves are NOT submitted here: the two shared solutions
#   permeability_formula/static_solve_out.e          (base mesh)
#   permeability_formula/static_solve_refined_out.e  (refined mesh)
# are expected to already exist / be uploaded before the elasticity jobs start.
#
# Use --with-static to also submit a per-case submit_static.sbatch first and
# chain via --dependency=afterok:<static_jid> (only the older self-contained
# cases have a per-case static script; the 5-case set shares the parent ones).
#
# Discovers every submit_elasticity.sbatch and submit_elasticity_refined.sbatch
# at depth <= 5 and submits each from its case directory so SLURM .o*/.e* logs
# and MOOSE output files land next to the input files.
#
# 5-case comparison set (each base + refined = 2 jobs), under permeability_formula/:
#   baseline/                control: damage_gradient normal, damage porosity
#   regularized_normal/      Group A: + regularize_crack_normal
#   strained-based/          Group A: principal_strain normal
#   porosity_bounded/        Group B: damage porosity, upper bound 0.065
#   porosity-strain-based/   Group B: porosity_update_model = strain
# Also discovered (base mesh only): the original top-level case, the
# undrained/drained EM sweeps, order_test/, viscosity/, etc. All share the two
# parent static_solve*.e files. Other 2d_hydromech subdirs (e.g. benchmark/)
# are NOT submitted by this script.
#
# Usage:
#   ./submit_all.sh                 # submit elasticity only (with confirmation)
#   ./submit_all.sh --yes           # submit without confirmation
#   ./submit_all.sh --dry-run       # show what would be submitted
#   ./submit_all.sh --new           # ONLY the 5-case comparison set x 2 meshes = 10 jobs
#   ./submit_all.sh --new --yes     # ... the 10 jobs, no confirmation prompt
#   ./submit_all.sh --with-static   # also submit static (chain via afterok)
#   ./submit_all.sh undrained       # only undrained EM sweep
#   ./submit_all.sh drained         # only drained EM sweep
#   ./submit_all.sh em_0p040        # only the EM=0.040 cases (both drained and undrained)
#   ./submit_all.sh order_test      # only the mixed-order test
#
#   --new restricts to: baseline, regularized_normal, strained-based,
#   porosity_bounded, porosity-strain-based (base + refined mesh of each).
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
NEW_ONLY=0
PATTERN=""

# The 5 case directories added for the permeability comparison set (each ships a
# base + refined mesh job => 10 jobs). --new restricts submission to exactly these.
NEW_CASES=(baseline regularized_normal strained-based porosity_bounded porosity-strain-based)

while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)       DRY_RUN=1 ;;
        --yes|-y)        SKIP_CONFIRM=1 ;;
        --with-static)   WITH_STATIC=1 ;;
        --new)           NEW_ONLY=1 ;;
        -h|--help)
            sed -n '2,53p' "$0"
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

# Discover every elasticity job: both the base-mesh submit_elasticity.sbatch and
# the refined-mesh submit_elasticity_refined.sbatch (the 5-case comparison set
# ships both; older sweep cases ship only the base). -maxdepth 5 handles
# ./<case>/, ./<group>/<case>/, ...
# Portable array fill for bash 3.2+ (mapfile is bash 4+ only).
ELAST_SCRIPTS=()
while IFS= read -r s; do
    ELAST_SCRIPTS+=("$s")
done < <(find "$SEARCH_ROOT" -maxdepth 5 -type f \
             \( -name 'submit_elasticity.sbatch' -o -name 'submit_elasticity_refined.sbatch' \) \
             | sort)

if [[ -n "$PATTERN" ]]; then
    FILTERED=()
    for s in "${ELAST_SCRIPTS[@]}"; do
        if [[ "$s" == *"$PATTERN"* ]]; then
            FILTERED+=("$s")
        fi
    done
    ELAST_SCRIPTS=("${FILTERED[@]}")
fi

# --new: keep only jobs whose case directory is one of the 5 NEW_CASES
# (matched on the case-dir basename, so base + refined of each are both kept,
# while the sweep/order/viscosity/top-level cases are dropped) => 10 jobs.
if [[ $NEW_ONLY -eq 1 ]]; then
    FILTERED=()
    for s in "${ELAST_SCRIPTS[@]}"; do
        case_base="$(basename "$(dirname "$s")")"
        for nc in "${NEW_CASES[@]}"; do
            if [[ "$case_base" == "$nc" ]]; then
                FILTERED+=("$s")
                break
            fi
        done
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

echo "Found ${#ELAST_SCRIPTS[@]} job(s):"
for s in "${ELAST_SCRIPTS[@]}"; do
    tag="base"; [[ "$s" == *_refined.sbatch ]] && tag="refined"
    printf '  %-10s %s\n' "[$tag]" "$(dirname "${s#$SCRIPT_DIR/}")"
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
    elast_name="$(basename "$elast_script")"     # submit_elasticity[_refined].sbatch
    # Stage label + display name distinguish the two meshes of a case.
    if [[ "$elast_name" == *_refined.sbatch ]]; then
        stage="elasticity_refined"; disp_name="$case_name (refined)"
    else
        stage="elasticity"; disp_name="$case_name"
    fi
    # The elasticity input this sbatch actually runs (parsed from its '-i' arg),
    # so the static-IC check below targets the matching base/refined solution.
    run_input="$(grep -oE '\-i +[^[:space:]]+\.i' "$elast_script" 2>/dev/null | head -1 | sed -E 's/^-i +//')"
    run_input="$(basename "${run_input:-}")"
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
        # Sanity check: resolve the static-solve output path referenced by the
        # specific input this sbatch runs (SolutionUserObject 'mesh = ...e').
        # Base inputs point at ../static_solve_out.e, refined inputs at
        # ../static_solve_refined_out.e. Warn only if the file isn't there.
        static_ref=""
        if [[ -n "$run_input" && -f "$run_input" ]]; then
            static_ref=$(grep -Eo 'mesh = [^[:space:]]+\.e' "$run_input" | head -1 | sed 's|^mesh = ||')
        fi
        if [[ -n "$static_ref" ]]; then
            if [[ ! -f "$static_ref" ]]; then
                echo "  [WARN] $disp_name: referenced static-solve output '$static_ref' not found (generate/upload it before the job starts)."
            fi
        elif [[ ! -f "static_solve_out.e" ]]; then
            # No SolutionUserObject reference parsed; fall back to the
            # historical check against ./static_solve_out.e.
            echo "  [WARN] $disp_name: static_solve_out.e not found; elasticity job will fail until you upload it."
        fi
    fi

    if [[ -n "$static_jid" ]]; then
        OUT="$(sbatch --dependency=afterok:"$static_jid" "$elast_name" 2>&1)" || true
    else
        OUT="$(sbatch "$elast_name" 2>&1)" || true
    fi
    if [[ "$OUT" =~ Submitted\ batch\ job\ ([0-9]+) ]]; then
        JID="${BASH_REMATCH[1]}"
        dep_note="${static_jid:+(afterok:$static_jid)}"
        printf '  [OK]   %-50s %-18s jobid=%s  %s\n' "$disp_name" "$stage" "$JID" "$dep_note"
        SUBMITTED_IDS+=("$JID:$case_name:$stage")
    else
        printf '  [FAIL] %-50s %-18s sbatch: %s\n' "$disp_name" "$stage" "$OUT"
        FAILED+=("$case_name:$stage")
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
