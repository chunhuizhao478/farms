"""
MPI-accelerated post-processing of on-fault point time histories for multiple cases.

Run with: mpiexec -n <N> python post_onfault_mpi.py

Strategy
- Rank 0 reads the first snapshot to find closest row indices for each target point,
  determines column ordering, prepares/clears the output directory, and broadcasts
  indices/columns to all ranks.
- Time steps are partitioned across ranks (contiguous blocks).
- Each rank reads only its assigned CSV snapshots and collects rows for the
  pre-selected indices, then writes one partial CSV per target point into a
  per-rank temp folder.
- Rank 0 merges partials for each point into a single CSV (sorted by time) and
  cleans up temp folders.

This script mirrors the configuration of post_onfault.py but uses MPI for speed.
"""

from __future__ import annotations
import os
import shutil
import math
from typing import Dict, Tuple, List
import sys

import pandas as pd
import numpy as np
from tqdm import tqdm

try:
    from mpi4py import MPI
except Exception as e:  # pragma: no cover
    raise SystemExit("mpi4py is required. Install with 'pip install mpi4py' and run with mpiexec.")

###############################
# User inputs (copy/adapt from post_onfault.py)
###############################

# Shared target points (x, y, z) used for all cases
# Edit as needed
TARGET_POINTS: List[Tuple[float, float, float]] = [
    (-12000, 0, -500),
    (-8000, 0, -500),
    (-4000, 0, -500),
    (0, 0, -500),
    (4000, 0, -500),
    (8000, 0, -500),
    (12000, 0, -500),
    (16000, 0, -500), #
    (-12000, 0, -1000),
    (-8000, 0, -1000),
    (-4000, 0, -1000),
    (0, 0, -1000),
    (4000, 0, -1000),
    (8000, 0, -1000),
    (12000, 0, -1000), #
    (-12000, 0, -1500),
    (-8000, 0, -1500),
    (-4000, 0, -1500),
    (0, 0, -1500),
    (4000, 0, -1500),
    (8000, 0, -1500),
    (12000, 0, -1500), #
]

# Case list (copy/adapt from your post_onfault.py)
COMMON_PATH = "/scratch1/10024/zhaochun/projects/farms_benchmark_08282025/development/paper_draft_depth_dependent/case1g/dynamic_solve/"
CASES = [
    {
        "name": "case1g_elastic_csv_on_fault", # name: A label used for outputs and logs.
        "data_dir": COMMON_PATH,
        "file_prefix": "dynamic_solve_elastic_csv_on_fault_", # file_prefix: The literal prefix of the input CSV files to read.
        "index": {"start": 2, "end": 1960, "step": 2, "pad": 4},
        "dt": 0.005,
        "output_dir": COMMON_PATH + "../postprocess/case1g_elastic_csv_on_fault/",
    }
]

# Optional global output root. If set (non-None), outputs will be written to
# os.path.join(output_root, case_name). A per-case 'output_dir' overrides this.
OUTPUT_ROOT = None  # e.g., "./point_time_histories"
# Whether to clear existing files inside the output directory before each run.
# Can be overridden per case by setting case_cfg["clean_output"] = True/False.
CLEAN_OUTPUT_DEFAULT = False


def resolve_path(path_str: str | None) -> str | None:
    if path_str is None:
        return None
    cleaned = path_str.replace("\\ ", " ")
    return os.path.normpath(os.path.expanduser(cleaned))


def distance(p1, p2):
    return math.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2 + (p1[2] - p2[2]) ** 2)


def partition_steps(steps: List[int], rank: int, size: int) -> List[int]:
    """Contiguous partition of steps across ranks."""
    n = len(steps)
    if size <= 1 or n == 0:
        return steps
    base = n // size
    rem = n % size
    start = rank * base + min(rank, rem)
    count = base + (1 if rank < rem else 0)
    return steps[start : start + count]


def partition_list(items: List, rank: int, size: int) -> List:
    """Contiguous partition of a list across ranks."""
    n = len(items)
    if size <= 1 or n == 0:
        return items
    base = n // size
    rem = n % size
    start = rank * base + min(rank, rem)
    count = base + (1 if rank < rem else 0)
    return items[start : start + count]


def process_case_mpi(case_cfg: dict, target_points: List[Tuple[float, float, float]]):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    name = case_cfg["name"]
    data_dir = resolve_path(case_cfg["data_dir"]) or case_cfg["data_dir"]
    file_prefix = case_cfg["file_prefix"]
    index_cfg = case_cfg.get("index", {})
    start = int(index_cfg.get("start", 0))
    end = int(index_cfg.get("end", -1))
    step = int(index_cfg.get("step", 1))
    pad = int(index_cfg.get("pad", 4))
    dt = float(case_cfg.get("dt", 1.0))

    # Output dir selection
    case_output_dir = case_cfg.get("output_dir")
    if case_output_dir:
        output_dir = resolve_path(case_output_dir)
    elif OUTPUT_ROOT:
        output_dir = os.path.join(resolve_path(OUTPUT_ROOT), name)
    else:
        output_dir = os.path.join(data_dir, "point_time_histories", name)

    # User option: whether to clear existing files in output_dir
    clean_output = bool(case_cfg.get("clean_output", CLEAN_OUTPUT_DEFAULT))

    # Rank 0 prepares/clears output dir conditionally
    if rank == 0:
        if os.path.exists(output_dir):
            if not os.path.isdir(output_dir):
                if clean_output:
                    try:
                        os.unlink(output_dir)
                    except Exception:
                        pass
                    os.makedirs(output_dir, exist_ok=True)
                else:
                    print(
                        f"[ERROR] Output path exists and is not a directory: {output_dir}. "
                        f"Set clean_output=True or change output_dir."
                    )
                    return
            else:
                if clean_output:
                    for entry in tqdm(os.listdir(output_dir), desc=f"{name}: clean output", disable=(rank != 0), file=sys.stdout):
                        path = os.path.join(output_dir, entry)
                        try:
                            if os.path.isfile(path) or os.path.islink(path):
                                os.unlink(path)
                            elif os.path.isdir(path):
                                shutil.rmtree(path)
                        except Exception:
                            pass
                # else: keep existing files
        else:
            os.makedirs(output_dir, exist_ok=True)
    comm.Barrier()

    # Determine time steps
    time_steps = list(range(start, end + 1, step))
    if not time_steps:
        if rank == 0:
            print(f"[WARN] No time steps for case '{name}' (start={start}, end={end}, step={step}). Skipping.")
        return

    # Parallel nearest-point search (vectorized) and column determination
    if rank == 0:
        first_file = f"{file_prefix}{start:0{pad}d}.csv"
        first_file_path = os.path.join(data_dir, first_file)
        found_first = os.path.isfile(first_file_path)
    else:
        found_first = False

    found_first = comm.bcast(found_first, root=0)
    if not found_first:
        if rank == 0:
            print(f"[ERROR] First file not found for case '{name}': {os.path.join(data_dir, f'{file_prefix}{start:0{pad}d}.csv')}\n        Skipping case.")
        return

    if rank == 0:
        df_first = pd.read_csv(os.path.join(data_dir, f"{file_prefix}{start:0{pad}d}.csv"))
        if df_first.empty:
            print(f"[WARN] No rows in first file for case '{name}'. Skipping case.")
            columns = None
            coords = None
        else:
            columns = list(df_first.columns)
            # Extract coordinates as numpy array for fast distance computation
            try:
                coords = df_first[["x", "y", "z"]].to_numpy(dtype=float)
            except Exception:
                # If coordinate columns missing, abort
                print(f"[ERROR] Columns x,y,z not found in '{name}' first file. Skipping case.")
                columns = None
                coords = None
    else:
        columns = None
        coords = None

    columns = comm.bcast(columns, root=0)
    coords = comm.bcast(coords, root=0)

    if columns is None or coords is None or coords.size == 0:
        return

    # Distribute target points across ranks for nearest-index search
    local_points = partition_list(target_points, rank, size)
    local_map: Dict[Tuple[float, float, float], int] = {}
    # Progress bar for point-finding on each rank (only show on rank 0 to reduce clutter)
    for tp in tqdm(
        local_points,
        desc=f"{name}: find points r{rank}",
        position=rank,
        leave=True,
        disable=False,
        file=sys.stdout,
    ):
        p = np.array(tp, dtype=float)
        diffs = coords - p  # (N,3)
        # Use squared distances to avoid sqrt
        d2 = np.einsum("ij,ij->i", diffs, diffs)
        idx = int(np.argmin(d2))
        local_map[tp] = idx

    gathered = comm.gather(local_map, root=0)
    if rank == 0:
        closest_point_indices: Dict[Tuple[float, float, float], int] = {}
        for part in gathered:
            if part:
                closest_point_indices.update(part)
    else:
        closest_point_indices = None

    closest_point_indices = comm.bcast(closest_point_indices, root=0)
    if closest_point_indices is None:
        return

    # Partition steps
    local_steps = partition_steps(time_steps, rank, size)

    # Each rank writes partials under a fresh tmp_rank{rank}
    tmp_dir = os.path.join(output_dir, f"tmp_rank_{rank}")
    if os.path.isdir(tmp_dir):
        shutil.rmtree(tmp_dir, ignore_errors=True)
    os.makedirs(tmp_dir, exist_ok=True)

    # In-memory collectors per point
    point_data: Dict[Tuple[float, float, float], List[List[float]]] = {p: [] for p in target_points}

    # Process assigned steps with a progress bar (only show on rank 0 to avoid clutter)
    step_iter = tqdm(
        local_steps,
        desc=f"{name}: steps r{rank}/{size}",
        position=rank,
        leave=True,
        disable=False,
        file=sys.stdout,
    )
    for t_step in step_iter:
        file_name = f"{file_prefix}{t_step:0{pad}d}.csv"
        file_path = os.path.join(data_dir, file_name)
        if not os.path.isfile(file_path):
            if rank == 0:
                print(f"[WARN] Missing file: {file_path}")
            continue
        df = pd.read_csv(file_path)
        t_val = t_step * dt
        for target_point, idx in closest_point_indices.items():
            if 0 <= idx < len(df):
                row = df.iloc[idx]
                point_data[target_point].append([t_val] + row.tolist())

    # Write one partial per point for this rank (with progress)
    items_to_write = [(p, rows) for p, rows in point_data.items() if rows]
    for point, rows in tqdm(
        items_to_write,
        desc=f"{name}: write r{rank}",
        total=len(items_to_write),
        position=rank,
        leave=True,
        disable=False,
        file=sys.stdout,
    ):
        x, y, z = point
        df_part = pd.DataFrame(rows, columns=["time"] + columns)
        part_name = f"part_rank{rank}_point_{x:.2f}_{y:.2f}_{z:.2f}.csv"
        df_part.to_csv(os.path.join(tmp_dir, part_name), index=False)

    comm.Barrier()

    # Parallel merge partials per point: distribute points across ranks
    local_points = partition_list(target_points, rank, size)
    for point in tqdm(
        local_points,
        desc=f"{name}: merge r{rank}",
        total=len(local_points),
        position=rank,
        leave=True,
        disable=False,
        file=sys.stdout,
    ):
        x, y, z = point
        final_name = f"{name}_point_{x:.2f}_{y:.2f}_{z:.2f}.csv"
        partials = []
        for r in range(size):
            pdir = os.path.join(output_dir, f"tmp_rank_{r}")
            pname = f"part_rank{r}_point_{x:.2f}_{y:.2f}_{z:.2f}.csv"
            ppath = os.path.join(pdir, pname)
            if os.path.isfile(ppath):
                try:
                    partials.append(pd.read_csv(ppath))
                except Exception:
                    pass
        if not partials:
            continue
        df_final = pd.concat(partials, ignore_index=True)
        df_final.sort_values("time", inplace=True)
        df_final.to_csv(os.path.join(output_dir, final_name), index=False)

    comm.Barrier()
    if rank == 0:
        # Cleanup tmp dirs (with progress)
        for r in tqdm(range(size), desc=f"{name}: cleanup tmp", disable=False, file=sys.stdout):
                pdir = os.path.join(output_dir, f"tmp_rank_{r}")
                if os.path.isdir(pdir):
                    shutil.rmtree(pdir, ignore_errors=True)
        print(f"[INFO] Case '{name}' complete. Output: {output_dir}")


def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if not CASES:
        if rank == 0:
            print("[INFO] No cases configured. Please edit CASES list in this script.")
        return

    for case_cfg in CASES:
        process_case_mpi(case_cfg, TARGET_POINTS)

    if rank == 0:
        print("\nAll cases complete.")


if __name__ == "__main__":
    main()
