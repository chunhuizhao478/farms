#!/usr/bin/env python3

"""
Post-process T1 traction output on the main fault and create a time-history plot.

The script expects CSV snapshots produced by test_me_maintb.i in the same folder,
with file names matching: test_me_maintb_csv_main_fault_XXXX.csv
where XXXX is the step number (every 40 solver steps).

Outputs:
  - t1_time_history_main_fault.png : heatmap of shear traction (MPa) vs arclength and time
  - t1_final_profile_main_fault.png : line plot of the final snapshot with comparison curves
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable, Tuple

import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
import pandas as pd


DATA_DIR = Path(__file__).resolve().parent
CSV_PATTERN = "test_me_maintb_csv_main_fault_*.csv"
ORIGIN = np.array([0.0, 0.036688], dtype=float)  # metres
DT = 3.015e-8  # seconds per solver step
OUTPUT_HEATMAP = DATA_DIR / "t1_time_history_main_fault.png"
OUTPUT_FINAL = DATA_DIR / "t1_final_profile_main_fault.png"
OUTPUT_MOVIE = DATA_DIR / "t1_main_fault_time_history.mp4"


def _collect_files(pattern: str) -> list[Path]:
    files = sorted(
        DATA_DIR.glob(pattern),
        key=lambda p: int(re.search(r"_(\d+)\.csv$", p.name).group(1)),
    )
    if not files:
        raise FileNotFoundError(
            f"No CSV files matching pattern '{pattern}' were found in {DATA_DIR}"
        )
    return files


def _unit_direction(df: pd.DataFrame) -> np.ndarray:
    rel = df[["x", "y"]].to_numpy() - ORIGIN
    norms = np.linalg.norm(rel, axis=1)
    if np.allclose(norms, 0):
        raise ValueError(
            "Cannot derive direction vector: all nodes coincide with the origin."
        )
    dir_vec = rel[norms.argmax()]
    return dir_vec / np.linalg.norm(dir_vec)


def _project_arclength(df: pd.DataFrame, direction: np.ndarray) -> np.ndarray:
    rel = df[["x", "y"]].to_numpy() - ORIGIN
    return rel @ direction


def _reference_arclength(df: pd.DataFrame, direction: np.ndarray) -> np.ndarray:
    df = df.copy()
    df["s"] = _project_arclength(df, direction)
    ref = df.sort_values("s").drop_duplicates(subset="s", keep="first")
    return ref["s"].to_numpy()


def _load_series(
    file: Path, direction: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    df = pd.read_csv(file)
    df["s"] = _project_arclength(df, direction)
    df = df.sort_values("s").drop_duplicates(subset="s", keep="first")
    return df["s"].to_numpy(), df["T1_aux"].to_numpy(), df["T2_aux"].to_numpy()


def _compute_time(step: int) -> float:
    return step * DT


def main() -> None:
  files = _collect_files(CSV_PATTERN)
  first_df = pd.read_csv(files[0])
  direction = _unit_direction(first_df)
  reference_s = _reference_arclength(first_df, direction)

  n_steps = len(files)
  n_points = reference_s.size
  t1_history = np.full((n_steps, n_points), np.nan, dtype=float)
  t2_history = np.full((n_steps, n_points), np.nan, dtype=float)
  times = np.zeros(n_steps, dtype=float)

  for idx, file in enumerate(files):
    step = int(re.search(r"_(\d+)\.csv$", file.name).group(1))
    times[idx] = _compute_time(step)

    s_values, t1_values, t2_values = _load_series(file, direction)

    # Remove potential duplicated arclengths that survived drop_duplicates due to floating tolerance.
    mask = np.diff(s_values, prepend=s_values[0] - 1e-12) > 0
    s_values = s_values[mask]
    t1_values = t1_values[mask]
    t2_values = t2_values[mask]

    if s_values.size < 2:
      raise ValueError(f"Insufficient data points after deduplication in file {file}")

    t1_diff_interp = np.interp(
        reference_s, s_values, t1_values, left=t1_values[0], right=t1_values[-1]
    )
    t2_diff_interp = np.interp(
        reference_s, s_values, t2_values, left=t2_values[0], right=t2_values[-1]
    )

    t1_history[idx, :] = t1_diff_interp
    t2_history[idx, :] = t2_diff_interp

  t1_history_mpa = -(t1_history / 1e6)  # convert to MPa and flip sign
  t2_history_mpa = t2_history / 1e6

  arc_lengths = reference_s
  times_us = times * 1e6  # microseconds

  arc_grid, time_grid = np.meshgrid(arc_lengths, times_us)

  fig, ax = plt.subplots(figsize=(10, 6))
  heat = ax.pcolormesh(
      arc_grid,
      time_grid,
      t1_history_mpa,
      shading="auto",
      cmap="RdBu_r",
  )
  ax.set_xlabel("Arclength from nucleation point (m)")
  ax.set_ylabel("Time (µs)")
  cbar = fig.colorbar(heat, ax=ax, pad=0.01)
  cbar.set_label("-T1_aux (MPa)")
  ax.set_title("Main Fault Shear Traction History")
  fig.tight_layout()
  fig.savefig(OUTPUT_HEATMAP, dpi=300)
  plt.close(fig)

  # Animated line plot over time (three curves every frame)
  shear_history = t1_history_mpa
  strength_mu_s = -0.7 * t2_history_mpa
  strength_mu_d = -0.1 * t2_history_mpa

  fig, ax = plt.subplots(figsize=(9, 4))
  line_shear, = ax.plot([], [], lw=2, label="-T1_aux")
  line_mu_s, = ax.plot([], [], lw=2, linestyle="--", label="-0.7 × T2_aux")
  line_mu_d, = ax.plot([], [], lw=2, linestyle="--", label="-0.1 × T2_aux")
  time_text = ax.text(0.02, 0.90, "", transform=ax.transAxes)
  ax.set_xlim(arc_lengths.min(), arc_lengths.max())
  combined_min = np.nanmin(np.vstack((shear_history, strength_mu_s, strength_mu_d)))
  combined_max = np.nanmax(np.vstack((shear_history, strength_mu_s, strength_mu_d)))
  ax.set_ylim(combined_min, min(combined_max, 20.0))
  ax.set_xlabel("Arclength from nucleation point (m)")
  ax.set_ylabel("Traction (MPa)")
  ax.grid(True, linestyle="--", alpha=0.4)
  ax.set_title("Main Fault Tractions vs Arclength")
  ax.legend(loc="upper right")

  def init_anim():
    line_shear.set_data([], [])
    line_mu_s.set_data([], [])
    line_mu_d.set_data([], [])
    time_text.set_text("")
    return line_shear, line_mu_s, line_mu_d, time_text

  def update_anim(frame: int):
    line_shear.set_data(arc_lengths, shear_history[frame])
    line_mu_s.set_data(arc_lengths, strength_mu_s[frame])
    line_mu_d.set_data(arc_lengths, strength_mu_d[frame])
    time_text.set_text(f"t = {times_us[frame]:.2f} µs")
    return line_shear, line_mu_s, line_mu_d, time_text

  fps = 30
  interval_ms = 1000 / fps
  anim = animation.FuncAnimation(
      fig,
      update_anim,
      frames=n_steps,
      init_func=init_anim,
      blit=True,
      interval=interval_ms,
  )

  try:
    anim.save(OUTPUT_MOVIE, writer="ffmpeg", dpi=200)
    print(f"Saved animation to {OUTPUT_MOVIE}")
  except (RuntimeError, OSError) as exc:
    print("WARNING: Failed to write MP4 movie (missing ffmpeg?).")
    print(f"         {exc}")
  plt.close(fig)

  # Plot the final snapshot for quick inspection
  fig, ax = plt.subplots(figsize=(9, 4))
  final_t1 = shear_history[-1]
  final_mu_s = strength_mu_s[-1]
  final_mu_d = strength_mu_d[-1]
  ax.plot(arc_lengths, final_t1, label=f"-T1_aux (t = {times_us[-1]:.2f} µs)")
  ax.plot(arc_lengths, final_mu_s, linestyle="--", label="-0.7 × T2_aux")
  ax.plot(arc_lengths, final_mu_d, linestyle="--", label="-0.1 × T2_aux")
  ax.set_xlabel("Arclength from nucleation point (m)")
  ax.set_ylabel("Traction (MPa)")
  ax.grid(True, linestyle="--", alpha=0.4)
  ax.legend()
  ax.set_title("Final Snapshot of Main Fault Tractions")
  ymin = min(np.nanmin(np.vstack((final_t1, final_mu_s, final_mu_d))), ax.get_ylim()[0])
  ax.set_ylim(ymin, 20.0)
  fig.tight_layout()
  fig.savefig(OUTPUT_FINAL, dpi=300)
  plt.close(fig)

  print(f"Saved heatmap to {OUTPUT_HEATMAP}")
  print(f"Saved final profile to {OUTPUT_FINAL}")


if __name__ == "__main__":
    main()
