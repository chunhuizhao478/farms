#!/bin/bash
# Bash script to generate 100 folders under the specified path

# Base directory where the folders will be created
BASE_DIR="/scratch1/10024/zhaochun/projects/farms_benchmark/examples/benchmark_tpv142D/fractal_stress"

# Loop to create folders output_1, output_2, ..., output_100
for i in {1..100}; do
    # Construct the folder name
    FOLDER="${BASE_DIR}/output_${i}"
    # Create the folder (including parent directories if they don't exist)
    mkdir -p "$FOLDER"
    echo "Created folder: $FOLDER"
done

echo "All folders have been created successfully!"