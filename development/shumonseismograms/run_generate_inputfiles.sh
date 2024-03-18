#!/bin/bash

# Path to your file
directory="."

# Directory to save the modified files
output_directory="./inputfiles"

# Ensure the output directory exists
mkdir -p "$output_directory"

echo "Starting file modifications..."

for i in {1..18}; do
    input_file_path="$directory/explicitdynamic.i"
    output_file_path="$output_directory/explicitdynamic_${i}.i"
    
    echo "Processing: $input_file_path -> $output_file_path"
    
    # Ensure the output file exists
    touch "$output_file_path"
    
    # Use sed to dynamically replace placeholders with the current iteration number
    # Note the escaping of special characters to match them literally in the search pattern
    sed "s/events_data_vx_num1/events_data_vx_num${i}/g" "$input_file_path" > "$output_file_path"
    sed -i '' "s/events_data_vy_num1/events_data_vy_num${i}/g" "$output_file_path"
    sed -i '' "s/events_data_time_num1/events_data_time_num${i}/g" "$output_file_path"
    sed -i '' "s|data_file = ./events_data_piecewisemultilinear|data_file = ../events_data_piecewisemultilinear|g" "$output_file_path"
    sed -i '' "s|files = ./events_data_piecewisemultilinear|files = ../events_data_piecewisemultilinear|g" "$output_file_path"

    echo "Updated file paths in $output_file_path for iteration $i"
done

echo "File modifications completed."