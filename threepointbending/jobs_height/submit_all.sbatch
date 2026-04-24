#!/bin/bash
# Submit all height parametric study jobs

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

for job in "$SCRIPT_DIR"/job_x*mm_y3.6mm.sbatch; do
  echo "Submitting $(basename "$job")"
  sbatch "$job"
done
