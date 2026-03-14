#!/bin/bash
# Submit all recover jobs for angle parametric study (x >= 1.50mm)

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

for job in "$SCRIPT_DIR"/job_x*mm_y1.6mm_recover.sbatch; do
  echo "Submitting $(basename "$job")"
  sbatch "$job"
done
