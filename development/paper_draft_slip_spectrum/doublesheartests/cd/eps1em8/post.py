#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

def guess_time_column(cols):
    lc = [c.lower() for c in cols]
    for k in ["time", "t", "Time"]:
        for c in cols:
            if c.lower() == k:
                return c
    # fallback: first column
    return cols[0]

def main():
    p = argparse.ArgumentParser(description="Plot columns from a MOOSE CSV output")
    p.add_argument("csv", help="CSV file (e.g. dynamic_solve_main_csv.csv)")
    p.add_argument("-c", "--columns", nargs="+",
                   help="Columns to plot (default: all except time)")
    p.add_argument("-o", "--out", help="Output image file (omit to show interactively)")
    p.add_argument("--rolling", type=int, default=0,
                   help="Apply rolling mean window (points)")
    p.add_argument("--ylim", nargs=2, type=float, help="y-axis limits")
    p.add_argument("--logy", action="store_true", help="Log y scale")
    p.add_argument("--separate", action="store_true",
                   help="Plot each variable in its own subplot")
    args = p.parse_args()

    path = Path(args.csv)
    if not path.is_file():
        print(f"File not found: {path}", file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(path)
    if df.empty:
        print("Empty CSV.", file=sys.stderr)
        sys.exit(1)

    time_col = guess_time_column(df.columns)
    if time_col not in df.columns:
        print("Time column not found.", file=sys.stderr)
        sys.exit(1)

    if args.columns:
        missing = [c for c in args.columns if c not in df.columns]
        if missing:
            print(f"Missing columns: {missing}", file=sys.stderr)
            sys.exit(1)
        plot_cols = args.columns
    else:
        plot_cols = [c for c in df.columns if c != time_col]

    t = df[time_col]

    if args.rolling > 1:
        df[plot_cols] = df[plot_cols].rolling(args.rolling, min_periods=1, center=True).mean()

    if args.separate:
        n = len(plot_cols)
        fig, axes = plt.subplots(n, 1, figsize=(8, 2.2 * n), sharex=True)
        if n == 1:
            axes = [axes]
        for ax, col in zip(axes, plot_cols):
            ax.plot(t, df[col], label=col)
            ax.set_ylabel(col)
            if args.logy:
                ax.set_yscale("log")
            if args.ylim:
                ax.set_ylim(args.ylim)
            ax.grid(True, alpha=0.3)
        axes[-1].set_xlabel(time_col)
    else:
        plt.figure(figsize=(9, 5))
        for col in plot_cols:
            plt.plot(t, df[col], label=col)
        if args.logy:
            plt.yscale("log")
        if args.ylim:
            plt.ylim(args.ylim)
        plt.xlabel(time_col)
        plt.ylabel("Value")
        plt.legend(fontsize="small", ncol=2)
        plt.grid(alpha=0.3)

    plt.tight_layout()
    if args.out:
        plt.savefig(args.out, dpi=200)
    else:
        plt.show()

if __name__ == "__main__":
    main()