#!/usr/bin/env bash
# Generate gold CSVs for the damaged_porosity material tests.
#
# Requires a working farms-opt executable (the conda moose env on macOS 26.x links the
# app dylib but NOT the final executable; run this where the executable links, e.g. CI /
# Linux). After running, verify the values against EXPECTED_VALUES.md before committing.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
APP="${FARMS_APP:-$HERE/../../../../farms-opt}"
mkdir -p "$HERE/gold"

run() {  # run <input> <produced_csv> <gold_name> [extra cli args...]
  local input="$1"; local produced="$2"; local gold="$3"; shift 3
  echo "==> $input"
  "$APP" -i "$HERE/$input" "$@"
  cp "$HERE/$produced" "$HERE/gold/$gold"
}

run damage_model.i  damage_model_out.csv        damage_model_out.csv
run strain_model.i  strain_model_out.csv        strain_model_out.csv
run strain_model.i  strain_lower_clamp_out.csv  strain_lower_clamp_out.csv \
    Materials/porosity_strain/porosity_lower_bound=0.05 Outputs/file_base=strain_lower_clamp_out

echo "Gold written to $HERE/gold. Confirm against EXPECTED_VALUES.md, then commit."
