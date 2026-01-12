import numpy as np
import os
import shutil

# process_line_sorted.py
# This script processes three input files:
# - res_ycoord_crack_10.txt (y-coordinate values)
# - res_xcoord_crack_10.txt (x-coordinate values)
# - res_ind_elem_crack_10.txt (element ids)
#
# It checks each line for points that lie within a tolerance of the line y = k*x + b.
# If exactly two points satisfy the condition, it computes the average x and y,
# records the element id from the corresponding line, and then sorts all results by average x.
# The element IDs and the subdomain IDs (generated as 100+i) are written as a single row in their output files.

# Parameters for the line equation and tolerance.
# k = 0.0      # Slope (change as needed)
# b = 0.0      # Intercept (change as needed)
# tol = 0.1    # Tolerance for checking if a point is on the line
# linename = "mf"
# subdomainids = 100

k = -0.5778      # Slope (change as needed)
b = 3461.9      # Intercept (change as needed)
tol = 10   # Tolerance for checking if a point is on the line
linename = "bf"
subdomainids = 300

subfolder = "preprocess_results_" + linename

if os.path.exists(subfolder):
    # Remove all files and subfolders inside the "result" folder.
    for filename in os.listdir(subfolder):
        file_path = os.path.join(subfolder, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)  # Remove a file or a symbolic link.
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)  # Remove a directory and all its contents.
        except Exception as e:
            print(f"Failed to delete {file_path}. Reason: {e}")
else:
    # Create the subfolder if it does not exist.
    os.makedirs(subfolder)

# List to store results as tuples: (avg_x, avg_y, elem_id)
results = []

# Open and read the three input files.
with open("raw_results/res_ycoord_crack_"+str(linename)+".txt", "r") as fy, \
     open("raw_results/res_xcoord_crack_"+str(linename)+".txt", "r") as fx, \
     open("raw_results/res_ind_elem_crack_"+str(linename)+".txt", "r") as fid:

    y_lines = fy.readlines()
    x_lines = fx.readlines()
    id_lines = fid.readlines()

    # Iterate over corresponding lines.
    for y_line, x_line, id_line in zip(y_lines, x_lines, id_lines):
        try:
            # Convert y and x lines into lists of floats.
            y_values = list(map(float, y_line.strip().split()))
            x_values = list(map(float, x_line.strip().split()))
        except ValueError:
            # Skip lines that cannot be parsed.
            continue

        # Identify indices of points that lie on the line within the tolerance.
        indices = [i for i, (x_val, y_val) in enumerate(zip(x_values, y_values))
                   if abs(y_val - (k * x_val + b)) <= tol]

        # If exactly two points satisfy the condition, process them.
        if len(indices) == 2:
            avg_x = sum(x_values[i] for i in indices) / 2
            avg_y = sum(y_values[i] for i in indices) / 2
            elem_id = id_line.strip()
            results.append((avg_x, avg_y, elem_id))

# Sort results by the computed average x value (ascending order).
results.sort(key=lambda tup: tup[0])

# Write the sorted (avg_x, avg_y) pairs to result.txt.
with open(subfolder+"/result_"+str(linename)+".txt", "w") as fout:
    for avg_x, avg_y, _ in results:
        fout.write(f"{avg_x} {avg_y}\n")

# Write the element ids in a single row (space separated) to result_ids.txt.
with open(subfolder+"/result_ids_"+str(linename)+".txt", "w") as fout_ids:
    # Extract element ids from results and join them with a space.
    fout_ids.write(" ".join(elem_id for _, _, elem_id in results))

# Write the subdomain IDs in a single row (space separated) to result_subdomainid.txt.
with open(subfolder+"/result_subdomainid_"+str(linename)+".txt", "w") as fout_subdomainid:
    # Generate subdomain IDs as 100 + index for each result.
    subdomain_ids = [str(subdomainids) for i in range(len(results))]
    fout_subdomainid.write(" ".join(subdomain_ids))