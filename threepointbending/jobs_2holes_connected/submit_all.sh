#!/bin/bash
# Submit all (d) 2 holes connected by thin line (diverging crack path) jobs.
# Jobs run in this folder: their ../../farms-opt and output files are relative to it.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR" || exit 1

for job in "$SCRIPT_DIR"/job_*.sbatch; do
  echo "Submitting $(basename "$job")"
  sbatch "$job"
done
