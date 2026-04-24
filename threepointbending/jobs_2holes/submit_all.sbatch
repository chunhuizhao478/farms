#!/bin/bash
# Submit all two-holes parametric study jobs

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

for job in "$SCRIPT_DIR"/job_x*mm_y1.6mm_radi*.sbatch; do
  echo "Submitting $(basename "$job")"
  sbatch "$job"
done
