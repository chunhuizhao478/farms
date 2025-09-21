#!/usr/bin/env python3
"""
Plot time vs velocity magnitude from one or more CSV files.

Expected columns in each CSV (case-insensitive):
- time (or t)
- maxvelx / velx / vx
- maxvely / vely / vy
- maxvelz / velz / vz / velz0

Usage:
  python post_multi.py file1.csv file2.csv ...
  python post_multi.py "*.csv"              # quotes for shell globbing on Windows

If no paths are provided, falls back to "dynamic_solve_main5_csv.csv" in the current directory.

Outputs:
- Single file: saves "<csv_dir>/time_vs_velmag.png"
- Multiple files: saves "<common_dir>/time_vs_velmag_combined.png"
"""
import sys
from pathlib import Path
import math
import pandas as pd
import matplotlib.pyplot as plt


def find_col(df, candidates):
    cols_lower = {c.lower(): c for c in df.columns}
    for name in candidates:
        if name.lower() in cols_lower:
            return cols_lower[name.lower()]
    return None


def load_time_velmag(csv_path: Path):
    """Return (time_series, vel_mag_series) for one csv. Raises on missing cols."""
    df = pd.read_csv(csv_path)
    time_col = find_col(df, ["time", "t"])
    vx_col = find_col(df, ["maxvelx", "velx", "vx"])
    vy_col = find_col(df, ["maxvely", "vely", "vy"])
    vz_col = find_col(df, ["maxvelz", "velz", "vz", "velz0"])

    if time_col is None or vx_col is None or vy_col is None or vz_col is None:
        raise ValueError(
            f"{csv_path.name}: missing required columns.\n"
            f"Found: {list(df.columns)}\n"
            "Needed: time/t and vel components (x,y,z)."
        )

    vel_mag = (df[vx_col] ** 2 + df[vy_col] ** 2 + df[vz_col] ** 2) ** 0.5
    return df[time_col], vel_mag


def main():
    # Resolve input list (allow zero, one, or many files; allow globs passed by the shell)
    args = sys.argv[1:]
    paths = [Path(a) for a in args]

    if not paths:
        # Fall back to the single-file behavior
        default_path = Path("dynamic_solve_main5_csv.csv")
        if not default_path.exists():
            print("No input provided and default CSV not found: dynamic_solve_main5_csv.csv", file=sys.stderr)
            sys.exit(2)
        paths = [default_path]

    # Expand any directories into *.csv within (non-recursive)
    expanded = []
    for p in paths:
        if p.is_dir():
            expanded.extend(sorted(p.glob("*.csv")))
        else:
            expanded.append(p)
    paths = expanded

    # Check existence
    missing = [p for p in paths if not p.exists()]
    if missing:
        for m in missing:
            print(f"CSV not found: {m}", file=sys.stderr)
        sys.exit(2)

    # Load and plot
    plt.figure()
    loaded_any = False
    for path in paths:
        try:
            t, vmag = load_time_velmag(path)
        except Exception as e:
            print(f"Skipping {path}: {e}", file=sys.stderr)
            continue
        label = path.stem
        plt.plot(t, vmag, linewidth=1.5, label=label)
        loaded_any = True

    if not loaded_any:
        print("No valid CSVs to plot.", file=sys.stderr)
        sys.exit(2)

    plt.xlabel("Time")
    plt.ylabel("Velocity magnitude")
    plt.yscale("log")
    plt.title("Time vs Velocity Magnitude")

    # Add legend if multiple series
    if len(paths) > 1:
        plt.legend()

    plt.tight_layout()

    # Save
    if len(paths) == 1:
        out_png = paths[0].with_name("time_vs_velmag.png")
    else:
        # Save to a common parent directory
        try:
            common = Path(Path(paths[0]).parent)
            for p in paths[1:]:
                # Walk up until common parent matches
                while not str(Path(p).parent).startswith(str(common)):
                    common = common.parent
            out_png = common / "time_vs_velmag_combined.png"
        except Exception:
            out_png = Path.cwd() / "time_vs_velmag_combined.png"

    plt.savefig(out_png, dpi=200)

    try:
        plt.show()
    except Exception:
        pass

    print(f"Saved figure to: {out_png.resolve()}")


if __name__ == "__main__":
    main()
