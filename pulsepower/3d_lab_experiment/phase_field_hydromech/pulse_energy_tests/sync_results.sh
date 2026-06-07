#!/bin/bash
# =============================================================================
# Pull result files from Frontera into the external-drive archive tree.
#
# Run this on your LOCAL machine (Mac), not on Frontera.
#
# Pulls case_* subdirectories under
#   3d_lab_experiment/phase_field_hydromech/pulse_energy_tests/
# Edit SUBDIRS below to change the selection.
#
# Default destination:
#   /Volumes/One Touch/Research/PulsePowerFracturing/3d_lab_experiment/phase_field_hydromech/pulse_energy_tests
# Override with --dest <path> or the LOCAL_DEST environment variable.
#
# Usage:
#   ./sync_results.sh                                         # defaults
#   ./sync_results.sh --dry-run                               # preview only
#   ./sync_results.sh zhaochun@frontera.tacc.utexas.edu       # explicit host
#   ./sync_results.sh --dest /some/other/path
#   LOCAL_DEST=/some/other/path ./sync_results.sh
#
# If no host is given, $FRONTERA_HOST is used, falling back to the ssh alias
# 'frontera' (configure in ~/.ssh/config).
#
# What is transferred (whitelist):
#   *.e     — Exodus output (main results)
#   *.csv   — CSV postprocessor output
#
# What is NOT transferred:
#   elasticity_mesh2x_out_fracture0.e (intermediate fracture-field dump),
#   checkpoint/*, .jitcache/*, *.i, SLURM *.o / *.e<jobid> logs, everything else.
#
# If you need checkpoints (e.g. to restart a run from the local box), add
# --include='*.cpr' --include='*.cpa' before the --exclude='*' line.
# =============================================================================

set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REMOTE_ROOT="/scratch2/10024/zhaochun/projects/farms_cdms_04192026/pulsepower/3d_lab_experiment/phase_field_hydromech/pulse_energy_tests"
DEFAULT_DEST="/Volumes/One Touch/Research/PulsePowerFracturing/3d_lab_experiment/phase_field_hydromech/pulse_energy_tests"

# Subdirectories under REMOTE_ROOT to sync. Each becomes a separate rsync pass.
SUBDIRS=(
    "case_E5J"
    "case_E5J_engtrack"
)

# Subdirectories where only *.csv should be transferred (skip *.e).
CSV_ONLY_SUBDIRS=(
    "case_E5J_engtrack"
)

# --- parse args (host, --dry-run, --dest <path>, in any order) ---
REMOTE=""
DRY=""
DEST_ARG=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)   DRY="--dry-run"; shift ;;
        --dest)      DEST_ARG="$2"; shift 2 ;;
        --dest=*)    DEST_ARG="${1#--dest=}"; shift ;;
        -h|--help)   sed -n '2,36p' "$0"; exit 0 ;;
        *)           REMOTE="$1"; shift ;;
    esac
done
REMOTE="${REMOTE:-${FRONTERA_HOST:-frontera}}"
LOCAL_DEST="${DEST_ARG:-${LOCAL_DEST:-$DEFAULT_DEST}}"

# Refuse to run when the destination lives under /Volumes/<drive> but that
# drive isn't mounted.
if [[ "$LOCAL_DEST" == /Volumes/* ]]; then
    vol_name="${LOCAL_DEST#/Volumes/}"
    vol_name="${vol_name%%/*}"
    if [[ ! -d "/Volumes/$vol_name" ]]; then
        echo "ERROR: /Volumes/$vol_name is not mounted; mount the drive or pass --dest <path>." >&2
        exit 1
    fi
fi

mkdir -p "$LOCAL_DEST" || {
    echo "ERROR: cannot create destination $LOCAL_DEST (check drive is writable)." >&2
    exit 1
}

echo "Remote: $REMOTE:$REMOTE_ROOT"
echo "Local:  $LOCAL_DEST"
echo "Subdirs: ${SUBDIRS[*]}"
echo "Mode:   ${DRY:-live}"
echo

for sub in "${SUBDIRS[@]}"; do
    # Decide whether to include *.e for this case.
    include_exodus=1
    for csvonly in "${CSV_ONLY_SUBDIRS[@]}"; do
        if [[ "$sub" == "$csvonly" ]]; then
            include_exodus=0
            break
        fi
    done

    if (( include_exodus )); then
        echo "==> Syncing $sub (csv + exodus)"
        mkdir -p "$LOCAL_DEST/$sub"
        rsync -avz --progress $DRY \
            --prune-empty-dirs \
            --include='*/' \
            --exclude='elasticity_mesh2x_out_fracture0.e' \
            --include='*.e' \
            --include='*.csv' \
            --exclude='checkpoint/***' \
            --exclude='.jitcache/***' \
            --exclude='*' \
            "$REMOTE:$REMOTE_ROOT/$sub/" \
            "$LOCAL_DEST/$sub/"
    else
        echo "==> Syncing $sub (csv only)"
        mkdir -p "$LOCAL_DEST/$sub"
        rsync -avz --progress $DRY \
            --prune-empty-dirs \
            --include='*/' \
            --include='*.csv' \
            --exclude='checkpoint/***' \
            --exclude='.jitcache/***' \
            --exclude='*' \
            "$REMOTE:$REMOTE_ROOT/$sub/" \
            "$LOCAL_DEST/$sub/"
    fi
    echo
done

echo "Done."
