import os
import numpy as np
import pandas as pd

# --- User settings ---
# Time stepping parameters
start = 20          # starting time step (e.g., "0000.csv")
end = 600      # ending time step
interval = 20      # time step interval
faultname = "main_fault" #"branch_fault" , "main_fault"  # name of the fault 

input_path = "/Users/chunhuizhao/projects/farms_benchmark/postprocess/tpv2053d_tet/test_csv/"
output_path = "/Users/chunhuizhao/projects/farms_benchmark/postprocess/tpv2053d_tet/farms_data_elemental_qp/"

# List of points of interest (x,y,z) coordinates
points_of_interest = [(     0, 0, -7500),
                      ( -4500, 0, -7500), 
                      ( -7500, 0, -7500), 
                      (-12000, 0, -7500)] #main fault

# File naming pattern: files are named with a four-digit time step
def get_filename(t, input_path, faultname):
    return f"{input_path}tpv2053D_tet_csv_{faultname}_{t:04d}.csv"

# --- Function definitions ---

def find_closest_index(df, target):
    """
    Given a DataFrame with 'x' and 'y' columns and a target coordinate,
    returns the row index of the point in df that is closest to the target.
    """
    # Ensure 'x' and 'y' columns are numeric
    df['x'] = pd.to_numeric(df['x'], errors='coerce')
    df['y'] = pd.to_numeric(df['y'], errors='coerce')
    df['z'] = pd.to_numeric(df['z'], errors='coerce')
    
    dx = df['x'] - target[0]
    dy = df['y'] - target[1]
    dz = df['z'] - target[2]
    # Use squared distance for efficiency
    dist2 = dx**2 + dy**2 + dz**2

    # Check if the distance series is not empty and contains non-NaN values
    if dist2.dropna().empty:
        raise ValueError("No valid numeric entries found in 'x' and 'y' columns to compute distances.")

    return dist2.idxmin()

# --- Main processing ---

# Create the list of time steps and filenames
time_steps = list(range(start, end + 1, interval))
filenames = [get_filename(t, input_path, faultname) for t in time_steps]

# Dictionaries to store time histories for each quantity and each point.
# Keys will be tuples: (point, quantity), e.g., ((-2000, 0), 'traction_x')
data_series = {}

# Define the list of quantities to extract
quantities = ['local_shear_jump', 'local_shear_jump_rate', 
              'local_shear_traction']

# Use the first file to determine the row indices for each point of interest.
if not os.path.exists(filenames[0]):
    raise FileNotFoundError(f"File {filenames[0]} not found. Check your file paths and naming convention.")

# Read the first file (assuming comma-separated values)
df_first = pd.read_csv(filenames[0])
# Convert x and y to numeric types
df_first['x'] = pd.to_numeric(df_first['x'], errors='coerce')
df_first['y'] = pd.to_numeric(df_first['y'], errors='coerce')
df_first['z'] = pd.to_numeric(df_first['y'], errors='coerce')

# print("Columns in the first file:", df_first.columns.tolist())
# print(df_first[['x', 'y']].head())
# print(df_first['x'])

# Map each point of interest to its closest row index (assuming mesh coordinates remain constant)
point_indices = {}
for pt in points_of_interest:
    idx = find_closest_index(df_first, pt)
    point_indices[pt] = idx
    print(f"Point {pt} closest row index: {idx}")

# Initialize time history lists for each quantity and each point.
for pt in points_of_interest:
    for qty in quantities:
        data_series[(pt, qty)] = {'time': [], 'value': []}

# Loop over each time step file and extract the desired quantities
for t in time_steps:
    fname = get_filename(t, input_path, faultname)
    if not os.path.exists(fname):
        print(f"Warning: File {fname} does not exist. Skipping.")
        continue

    # Read the CSV file (comma-separated)
    df = pd.read_csv(fname)
    # Ensure x and y are numeric
    df['x'] = pd.to_numeric(df['x'], errors='coerce')
    df['y'] = pd.to_numeric(df['y'], errors='coerce')
    df['z'] = pd.to_numeric(df['z'], errors='coerce')
    
    for pt in points_of_interest:
        idx = point_indices[pt]
        # Append the current time and value for each quantity
        for qty in quantities:
            if qty not in df.columns:
                raise ValueError(f"Column {qty} not found in file {fname}. Check the header names.")
            value = df.loc[idx, qty]
            data_series[(pt, qty)]['time'].append(t)
            data_series[(pt, qty)]['value'].append(value)

# --- Save the time histories to files ---

def make_output_filename(quantity, point):
    """
    Generate an output filename given the quantity and point of interest.
    The strike value is computed from the absolute x-coordinate in km.
    """
    strike_km = point[0] / 1000.0
    dip_km = point[2] / 1000.0
    strike_str = f"{strike_km:.1f}"  # One decimal place
    dip_str = f"{dip_km:.1f}"  # One decimal place

    return f"{quantity}_strike{strike_str}_dip{dip_str}.txt"

# Save each time series into a separate text file.
for key, ts in data_series.items():
    pt, qty = key
    out_fname = make_output_filename(qty, pt)
    # Create a DataFrame from only the quantity values
    df_out = pd.DataFrame(ts['value'], columns=[qty])
    out_full_path = os.path.join(output_path, out_fname)
    # Save to CSV without header and index, using tab as the separator
    df_out.to_csv(out_full_path, index=False, header=False, sep='\t')
    print(f"Saved time history for {qty} at point {pt} to file: {out_full_path}")

print("All time histories saved successfully.")

