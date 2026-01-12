"""
Post-processing point time histories for multiple cases.
Created By Chunhui Zhao, Aug 17th, 2025

This script supports processing multiple datasets ("cases"). For each case, you
can specify the data directory, file prefix, index range (start/end/step and
zero-padding), and time step dt. A shared list of target points is used across
all cases. Outputs are written under a case-specific subfolder to avoid
overwriting.
"""

import os
import shutil
import pandas as pd
import numpy as np
from tqdm import tqdm

###############################
# User inputs
###############################

# Shared target points (x, y, z) used for all cases
target_points = [
    (-24000, 4000, 0),
    (-20000, 4000, 0),
    (-16000, 4000, 0),
    (-12000, 4000, 0),
    (-8000, 4000, 0),
    (-4000, 4000, 0),
    (0, 4000, 0),
    (4000, 4000, 0),
    (8000, 4000, 0),
    (12000, 4000, 0),
    (16000, 4000, 0),
    (20000, 4000, 0),
    (24000, 4000, 0),
    (-24000, 4000, -10000),
    (-20000, -4000, -10000),
    (-16000, -4000, -10000),
    (-12000, -4000, -10000),
    (-8000, -4000, -10000),
    (-4000, -4000, -10000),
    (0, -4000, -10000),
    (4000, -4000, -10000),
    (8000, -4000, -10000),
    (12000, -4000, -10000),
    (16000, -4000, -10000),
    (20000, -4000, -10000),
    (24000, -4000, -10000),
]

# List of cases to process. For each case specify:
# - name: label used to segregate outputs
# - data_dir: path to the folder containing CSV snapshots
# - file_prefix: prefix of CSV files before the zero-padded index
# - index: dict with start, end, step, and optional pad (default 4)
# - dt: time step size

#case 2e
COMMON_PATH = "/scratch1/10024/zhaochun/projects/farms_benchmark_08282025/development/paper_draft_depth_dependent/case1g/dynamic_solve/"
cases = [
    {
        "name": "case1g_elastic_csv_off_fault",
        "data_dir": COMMON_PATH,
        "file_prefix": "dynamic_solve_elastic_csv_off_fault_",
        "index": {"start": 2, "end": 1960, "step": 2, "pad": 4},
        "dt": 0.005,
        "output_dir": COMMON_PATH + "../postprocess/case1g_elastic_csv_off_fault/",
    },
]

# Calculate distance between two 3D points
def distance(p1, p2):
    return np.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2)


# Optional global output root. If set (non-None), outputs will be written to
# os.path.join(output_root, case_name). A per-case 'output_dir' overrides this.
output_root = None  # e.g., "./point_time_histories"


def resolve_path(path_str: str) -> str:
    """Resolve a user-provided path string for Python FS usage.

    - Expand '~'
    - Convert shell-escaped spaces ('\\ ') to actual spaces
    - Normalize the path
    """
    if path_str is None:
        return None
    cleaned = path_str.replace("\\ ", " ")
    return os.path.normpath(os.path.expanduser(cleaned))

def process_case(case_cfg):
    name = case_cfg["name"]
    data_dir = resolve_path(case_cfg["data_dir"])
    file_prefix = case_cfg["file_prefix"]
    index_cfg = case_cfg.get("index", {})
    start = int(index_cfg.get("start", 0))
    end = int(index_cfg.get("end", -1))
    step = int(index_cfg.get("step", 1))
    pad = int(index_cfg.get("pad", 4))
    dt = float(case_cfg.get("dt", 1.0))

    # Prepare output directory per case
    case_output_dir = case_cfg.get("output_dir")
    if case_output_dir:
        output_dir = resolve_path(case_output_dir)
    elif output_root:
        output_dir = os.path.join(resolve_path(output_root), name)
    else:
        output_dir = os.path.join(data_dir, "point_time_histories", name)
    # Ensure output directory exists and is empty
    if os.path.exists(output_dir):
        if not os.path.isdir(output_dir):
            # If a file exists with same name, remove it and create directory
            try:
                os.unlink(output_dir)
            except Exception as e:
                print(f"[WARN] Could not remove existing file at output path '{output_dir}': {e}")
            os.makedirs(output_dir, exist_ok=True)
        else:
            # Clear existing contents
            for entry in os.listdir(output_dir):
                path = os.path.join(output_dir, entry)
                try:
                    if os.path.isfile(path) or os.path.islink(path):
                        os.unlink(path)
                    elif os.path.isdir(path):
                        shutil.rmtree(path)
                except Exception as e:
                    print(f"[WARN] Failed to remove '{path}': {e}")
    else:
        os.makedirs(output_dir, exist_ok=True)

    # Time steps
    time_steps = list(range(start, end + 1, step))
    time_values = [i * dt for i in time_steps]

    if not time_steps:
        print(f"[WARN] No time steps for case '{name}' (start={start}, end={end}, step={step}). Skipping.")
        return

    # First file to anchor closest-point mapping
    first_file = f"{file_prefix}{start:0{pad}d}.csv"
    first_file_path = os.path.join(data_dir, first_file)

    if not os.path.isfile(first_file_path):
        print(f"[ERROR] First file not found for case '{name}': {first_file_path}")
        return

    print(f"\n=== Case: {name} ===")
    print("Finding closest points to targets...")
    df_first = pd.read_csv(first_file_path)

    # Map each target point to the closest row index in the first snapshot
    closest_point_indices = {}
    for target_point in target_points:
        distances = []
        for idx, row in df_first.iterrows():
            current_point = (row["x"], row["y"], row["z"])
            dist = distance(current_point, target_point)
            distances.append((idx, dist, current_point))

        if not distances:
            print(f"[WARN] No rows in first file for case '{name}'. Skipping case.")
            return

        closest_idx, min_dist, actual_point = min(distances, key=lambda x: x[1])
        closest_point_indices[target_point] = closest_idx

        print(f"Target: {target_point}")
        print(f"Closest: {actual_point}")
        print(f"Distance: {min_dist:.2f}")
        print("-" * 30)

    # Collect data per target point
    point_data = {point: [] for point in target_points}
    columns = None

    print("Processing time steps...")
    for i, t_step in tqdm(
        list(enumerate(time_steps)), total=len(time_steps), desc=f"{name}: time steps"
    ):
        file_name = f"{file_prefix}{t_step:0{pad}d}.csv"
        file_path = os.path.join(data_dir, file_name)

        if not os.path.isfile(file_path):
            tqdm.write(f"Missing file: {file_path}")
            continue

        df = pd.read_csv(file_path)
        if columns is None:
            columns = list(df.columns)

        for target_point, idx in closest_point_indices.items():
            if 0 <= idx < len(df):
                row = df.iloc[idx]
                point_data[target_point].append([time_values[i]] + row.tolist())

    # Write outputs
    print("Writing output files...")
    for point in tqdm(point_data.keys(), desc=f"{name}: writing files"):
        data = point_data[point]
        if not data:
            tqdm.write(f"No data found for point {point} (case '{name}')")
            continue

        x, y, z = point
        df_point = pd.DataFrame(data, columns=["time"] + (columns or []))
        # Include case name in filename to be explicit, even within subfolder
        filename = f"{name}_point_{x:.2f}_{y:.2f}_{z:.2f}.csv"
        df_point.to_csv(os.path.join(output_dir, filename), index=False)
        tqdm.write(f"[{name}] Time history for point near ({x}, {y}, {z}) written.")


def main():
    if not cases:
        print("[INFO] No cases configured. Please edit 'cases' list in this script.")
        return
    for case_cfg in cases:
        process_case(case_cfg)
    print("\nAll cases complete.")


if __name__ == "__main__":
    main()