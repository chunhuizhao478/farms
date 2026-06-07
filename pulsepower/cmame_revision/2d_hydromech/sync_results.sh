#!/bin/bash
# =============================================================================
# Pull result files from Frontera into the external-drive archive tree.
#
# Run this on your LOCAL machine (Mac), not on Frontera.
#
# Only the hydromech case subdirectories are pulled (edit SUBDIRS below to
# change the selection; default is every case currently in this folder).
#
# Default destination:
#   /Volumes/One Touch/Research/PulsePowerFracturing/cmame_revision/2d_hydromech
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
#   *.e     — Exodus output (main results, including static_solve_out.e)
#   *.csv   — CSV postprocessor output
#
# What is NOT transferred:
#   elasticity_out_fracture0.e (intermediate fracture-subapp dump),
#   checkpoint/*, .jitcache/*, *.i, SLURM *.o / *.e<jobid> logs, everything else.
#
# If you need checkpoints (e.g. to restart a run from the local box), add
# --include='*.cpr' --include='*.cpa' before the --exclude='*' line.
# =============================================================================

set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REMOTE_ROOT="/scratch2/10024/zhaochun/projects/farms_cdms_04192026/pulsepower/cmame_revision/2d_hydromech"
DEFAULT_DEST="/Volumes/One Touch/Research/PulsePowerFracturing/cmame_revision/2d_hydromech"

# Subdirectories under REMOTE_ROOT to sync. Each becomes a separate rsync pass.
# Currently set to the 5-case permeability-formula comparison set (each case has
# a baseline elasticity_E1d25.i and a refined elasticity_E1d25_refined.i twin).
SUBDIRS=(
    "permeability_formula/baseline"
    "permeability_formula/porosity-strain-based"
    "permeability_formula/porosity_bounded"
    "permeability_formula/regularized_normal"
    "permeability_formula/strained-based"
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
# drive isn't mounted. Otherwise mkdir -p silently creates stub directories
# on the internal disk under /Volumes and rsync then fails with a confusing
# "Permission denied" when SIP blocks writes to the nested path.
if [[ "$LOCAL_DEST" == /Volumes/* ]]; then
    vol_name="${LOCAL_DEST#/Volumes/}"
    vol_name="${vol_name%%/*}"
    if [[ ! -d "/Volumes/$vol_name" ]]; then
        echo "ERROR: /Volumes/$vol_name is not mounted; mount the drive or pass --dest <path>." >&2
        exit 1
    fi
fi

# Make sure the destination exists; rsync will create subdirs but not the root.
mkdir -p "$LOCAL_DEST" || {
    echo "ERROR: cannot create destination $LOCAL_DEST (check drive is writable)." >&2
    exit 1
}

echo "Remote: $REMOTE:$REMOTE_ROOT"
echo "Local:  $LOCAL_DEST"
echo "Subdirs: ${SUBDIRS[*]}"
echo "Mode:   ${DRY:-live}"
echo

# rsync rules: allow directory traversal, whitelist Exodus + CSV, reject rest.
# Order matters: the first matching include/exclude wins, and --exclude='*'
# at the end drops anything not already whitelisted.
for sub in "${SUBDIRS[@]}"; do
    echo "==> Syncing $sub"
    mkdir -p "$LOCAL_DEST/$sub"
    rsync -avz --progress $DRY \
        --prune-empty-dirs \
        --include='*/' \
        --exclude='elasticity_out_fracture0.e' \
        --include='*.e' \
        --include='*.csv' \
        --exclude='checkpoint/***' \
        --exclude='.jitcache/***' \
        --exclude='*' \
        "$REMOTE:$REMOTE_ROOT/$sub/" \
        "$LOCAL_DEST/$sub/"
    echo
done

echo "Done."
