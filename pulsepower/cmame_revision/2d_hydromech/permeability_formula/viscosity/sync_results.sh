#!/bin/bash
# =============================================================================
# Pull viscosity/ result files from Frontera into the external-drive archive.
#
# Run this on your LOCAL machine (Mac), not on Frontera.
#
# Default destination:
#   /Volumes/One Touch/Research/PulsePowerFracturing/cmame_revision/2d_hydromech/permeability_formula/viscosity
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
# Cases pulled:
#   pure_solid/   mu_1em4/   mu_1em3/   mu_1em1/
#
# What is transferred (whitelist):
#   *.e     — Exodus output
#   *.csv   — CSV postprocessor output
#
# What is NOT transferred:
#   *_fracture0.e (intermediate fracture-field dump),
#   *_checkpoint_cp/* (checkpoint dirs), .jitcache/*, *.i, SLURM logs.
# =============================================================================

set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REMOTE_ROOT="/scratch2/10024/zhaochun/projects/farms_cdms_04192026/pulsepower/cmame_revision/2d_hydromech/permeability_formula/viscosity"
DEFAULT_DEST="/Volumes/One Touch/Research/PulsePowerFracturing/cmame_revision/2d_hydromech/permeability_formula/viscosity"

# Subdirectories under REMOTE_ROOT to sync. Each becomes a separate rsync pass.
SUBDIRS=(
    "pure_solid"
    "mu_1em4"
    "mu_1em3"
    "mu_1em1"
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
        -h|--help)   sed -n '2,30p' "$0"; exit 0 ;;
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
    echo "==> Syncing $sub"
    mkdir -p "$LOCAL_DEST/$sub"
    rsync -avz --progress $DRY \
        --prune-empty-dirs \
        --include='*/' \
        --exclude='*_fracture0.e' \
        --exclude='*_checkpoint_cp/***' \
        --exclude='.jitcache/***' \
        --include='*.e' \
        --include='*.csv' \
        --exclude='*' \
        "$REMOTE:$REMOTE_ROOT/$sub/" \
        "$LOCAL_DEST/$sub/"
    echo
done

echo "Done."
