#!/bin/bash
# =============================================================================
# Pull result files from Frontera back into the local 2d_solid/ tree.
#
# Run this on your LOCAL machine, not on Frontera.
#
# Usage:
#   ./sync_results.sh                              # use $FRONTERA_HOST or 'frontera'
#   ./sync_results.sh chunhui@frontera.tacc.utexas.edu
#   ./sync_results.sh --dry-run                    # show what would transfer
#   ./sync_results.sh <host> --dry-run             # both
#
# Assumes you have an ssh entry for Frontera (e.g. in ~/.ssh/config) named
# either $FRONTERA_HOST or "frontera" and the scratch tree layout matches
# the absolute paths hard-coded in the SBATCH scripts.
#
# What is transferred (whitelist):
#   *.e                   — Exodus output (main results)
#   *.csv                 — CSV postprocessor output
#   submitted_jobs.txt    — the submission log written by submit_all.sh
#   *.o[0-9]* *.e[0-9]*   — SLURM stdout / stderr log files
#
# What is NOT transferred (so the pull stays small and fast):
#   checkpoint/*          — MOOSE checkpoint files (often multi-GB)
#   .jitcache/            — MOOSE AD JIT caches
#   *.i                   — input files (they live in git already)
#
# If you need checkpoints (e.g. to restart a run from the local box), add
# --include='*.cpr' --include='*.cpa' before the --exclude='*' line.
# =============================================================================

set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REMOTE_ROOT="/scratch2/10024/zhaochun/projects/farms_cdms_04192026/pulsepower/cmame_revision/2d_solid"

# --- parse args (host and/or --dry-run, in any order) ---
REMOTE=""
DRY=""
for arg in "$@"; do
    case "$arg" in
        --dry-run) DRY="--dry-run" ;;
        -h|--help) sed -n '2,32p' "$0"; exit 0 ;;
        *)         REMOTE="$arg" ;;
    esac
done
REMOTE="${REMOTE:-${FRONTERA_HOST:-frontera}}"

echo "Remote: $REMOTE:$REMOTE_ROOT"
echo "Local:  $SCRIPT_DIR"
echo "Mode:   ${DRY:-live}"
echo

# rsync rules: allow directory traversal, whitelist result files, reject rest.
# Order matters: the first matching include/exclude wins, and --exclude='*'
# at the end drops anything not already whitelisted.
rsync -avz --progress $DRY \
    --prune-empty-dirs \
    --include='*/' \
    --include='*.e' \
    --include='*.csv' \
    --include='submitted_jobs.txt' \
    --include='*.o[0-9]*' \
    --include='*.e[0-9]*' \
    --exclude='checkpoint/***' \
    --exclude='.jitcache/***' \
    --exclude='*' \
    "$REMOTE:$REMOTE_ROOT/" \
    "$SCRIPT_DIR/"

echo
echo "Done."
