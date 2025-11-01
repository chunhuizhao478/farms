#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Scan CSV files for two conditions (Cond1: closest node jump‐rate > threshold; "
            "Cond2: any node with x ≥ x_thresh has jump‐rate > threshold). "
            "When a condition is met, save separate plots for local_jump and local_jump_rate "
            "(each with a vertical line at the flip_index arc‐length) and export CSV data "
            "containing arc_length vs. local_jump and local_jump_rate. "
            "Finally, save a summary CSV listing the first step each condition was triggered."
        )
    )

    # 1) Folder / filename parameters
    parser.add_argument(
        "-p", "--P_cases",
        nargs="+",
        default=["case_P15"],
        help="List of P_case folders (e.g. case_P15 case_P20)."
    )
    parser.add_argument(
        "-d", "--delta_cases",
        nargs="+",
        default=["case_deltasigma3"],
        help="List of delta_case folders (e.g. case_deltasigma3 case_deltasigma4)."
    )
    parser.add_argument(
        "-s", "--stem",
        type=str,
        default="test_me_maintb_csv_main_fault_",
        help="Filename stem before the step number (default: 'test_me_maintb_csv_main_fault_')."
    )

    # 2) Geometry & threshold parameters
    parser.add_argument(
        "--x0",
        type=float,
        default=0.140641,
        help="X coordinate of the target point for Condition 1 (default: 0.140641)."
    )
    parser.add_argument(
        "--y0",
        type=float,
        default=0.115434,
        help="Y coordinate of the target point for Condition 1 (default: 0.115434)."
    )
    parser.add_argument(
        "--x_thresh",
        type=float,
        default=0.142127,
        help="X threshold for Condition 2 (default: 0.142127)."
    )
    parser.add_argument(
        "--rate_thresh",
        type=float,
        default=3.0,
        help="Jump‐rate threshold (m/s) for both conditions (default: 3.0)."
    )

    # 3) Flip index + rotation
    parser.add_argument(
        "--flip_index",
        type=int,
        default=336,
        help="Row index from which to flip signs of jump columns (default: 336)."
    )
    parser.add_argument(
        "--rot_deg",
        type=float,
        default=29.0,
        help="Rotation angle in degrees for projecting jump onto tangent (default: 29°)."
    )

    # 4) Step scanning range
    parser.add_argument(
        "--start_step",
        type=int,
        default=40,
        help="First step number to check (default: 40)."
    )
    parser.add_argument(
        "--end_step",
        type=int,
        default=16000,
        help="Last step number to check (default: 16000)."
    )
    parser.add_argument(
        "--step_increment",
        type=int,
        default=40,
        help="Increment between steps (default: 40)."
    )

    # 5) Output and verbosity
    parser.add_argument(
        "-o", "--output_csv",
        type=str,
        default="summary.csv",
        help="Name of the output CSV summary (default: summary.csv)."
    )
    parser.add_argument(
        "--plots_dir",
        type=str,
        default="plots",
        help="Directory in which to save plot PNGs (default: ./plots)."
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Print detailed progress messages."
    )

    return parser.parse_args()


def ensure_directory(path):
    if not os.path.isdir(path):
        os.makedirs(path, exist_ok=True)


def plot_jump_and_save(
    arclen: np.ndarray,
    local_jump_or_rate: np.ndarray,
    highlight_idx: int,
    step: int,
    P: str,
    delta: str,
    condition_label: str,
    flip_arc: float,
    y_label: str,
    curve_label: str,
    plots_dir: str
):
    """
    Create and save a single plot of arc length vs. local_jump_or_rate,
    highlight the point where the condition was triggered, draw a vertical line at flip_arc,
    and save as a PNG file under plots_dir. Returns the saved PNG path.
    """
    plt.figure(figsize=(10, 5))
    plt.plot(arclen, local_jump_or_rate, marker='o', linestyle='-', label=curve_label)

    # Draw vertical line at flip_arc
    if flip_arc is not None:
        plt.axvline(flip_arc, color='gray', linestyle='--', label=f'flip_index @ {flip_arc:.3f} m')

    # Highlight the triggering point
    xh = arclen[highlight_idx]
    yh = local_jump_or_rate[highlight_idx]
    plt.scatter(xh, yh, s=100, c='k', zorder=5)

    # Annotate with arclen and y-value
    label_text = f"(arclen={xh:.3f}, {y_label}={yh:.2f})"
    plt.annotate(
        label_text,
        xy=(xh, yh),
        xytext=(xh + 0.01, yh + np.sign(yh)*0.2),
        arrowprops=dict(arrowstyle="->")
    )

    plt.xlabel("Arc length (m)")
    plt.ylabel(y_label)
    plt.title(f"{P}/{delta} — {condition_label} @ step {step:04d}\n{curve_label} triggered")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    # Choose filename based on whether this is jump or jump_rate
    suffix = "jump" if y_label.lower().startswith("jump ") else "jump_rate"
    filename = f"{P}_{delta}_{condition_label}_step{step:04d}_{suffix}.png"
    save_path = os.path.join(plots_dir, filename)
    plt.savefig(save_path, dpi=150)
    plt.close()

    return save_path


def scan_case(
    root_dir: str,
    P: str,
    delta: str,
    stem: str,
    x0: float,
    y0: float,
    x_thresh: float,
    rate_thresh: float,
    flip_index: int,
    rot_rad: float,
    start_step: int,
    end_step: int,
    step_inc: int,
    plots_dir: str,
    data_dir: str,
    verbose: bool
) -> (int, int, list, list, list, list):
    """
    Scan all CSVs for a given (P, delta). Return:
      (cond1_step, cond2_step,
       saved_jump_plots, saved_jump_rate_plots,
       saved_data_files, saved_rate_data_files)
    where each step is the first at which the respective condition was met (or None),
    saved_jump_plots is a list of PNG paths for the local_jump plots,
    saved_jump_rate_plots is for local_jump_rate plots,
    and saved_data_files is a list of CSV paths for data 
    (arclen, local_jump) & saved_rate_data_files is a list of CSV paths for (arclen, local_jump_rate).
    """
    cond1_step = None
    cond2_step = None
    saved_jump_plots = []
    saved_jump_rate_plots = []
    saved_jump_data = []
    saved_rate_data = []

    folder = os.path.join(root_dir, P, delta)
    if verbose:
        print(f"\n=== Scanning folder: {folder} ===")

    for step in range(start_step, end_step + 1, step_inc):
        fname = f"{stem}{step:04d}.csv"
        fpath = os.path.join(folder, fname)
        if not os.path.isfile(fpath):
            continue

        df = pd.read_csv(fpath)

        # -------------------------
        # 1) Flip sign for rows ≥ flip_index
        # -------------------------
        df.loc[flip_index:, 'jump_x_aux']      *= -1
        df.loc[flip_index:, 'jump_x_rate_aux'] *= -1
        df.loc[flip_index:, 'jump_y_aux']      *= -1
        df.loc[flip_index:, 'jump_y_rate_aux'] *= -1

        # -------------------------
        # 2) Compute arc length
        # -------------------------
        dx = df['x'] - df['x'].iloc[0]
        dy = df['y'] - df['y'].iloc[0]
        arclen = np.sqrt(dx**2 + dy**2)

        # Determine arc‐length at flip_index (if valid)
        flip_arc = arclen[flip_index] if flip_index < len(arclen) else None

        # -------------------------
        # 3) Compute local_jump_rate in the tangent direction
        # -------------------------
        local_jump_rate = (
            df['jump_x_rate_aux'] * np.cos(rot_rad)
          + df['jump_y_rate_aux'] * np.sin(rot_rad)
        )

        # -------------------------
        # 4) Compute local_jump (not rate) in the tangent direction
        # -------------------------
        local_jump = (
            df['jump_x_aux'] * np.cos(rot_rad)
          + df['jump_y_aux'] * np.sin(rot_rad)
        )

        # -----------------------------------
        # 5) Check Condition 1: closest node
        # -----------------------------------
        if cond1_step is None:
            dist = np.hypot(df['x'] - x0, df['y'] - y0)
            idx0 = dist.idxmin()
            rate0 = local_jump_rate.iloc[idx0]
            if rate0 > rate_thresh:
                cond1_step = step
                if verbose:
                    print(f"[{P}/{delta}] Cond1 met @ step {step:04d}, rate = {rate0:.2f} m/s")
                # Plot local_jump
                pj = plot_jump_and_save(
                    arclen=arclen,
                    local_jump_or_rate=local_jump,
                    highlight_idx=idx0,
                    step=step,
                    P=P,
                    delta=delta,
                    condition_label="Cond1",
                    flip_arc=flip_arc,
                    y_label="jump (m)",
                    curve_label="tangent jump",
                    plots_dir=plots_dir
                )
                saved_jump_plots.append(pj)
                # Plot local_jump_rate
                pr = plot_jump_and_save(
                    arclen=arclen,
                    local_jump_or_rate=local_jump_rate,
                    highlight_idx=idx0,
                    step=step,
                    P=P,
                    delta=delta,
                    condition_label="Cond1",
                    flip_arc=flip_arc,
                    y_label="jump_rate (m/s)",
                    curve_label="tangent jump rate",
                    plots_dir=plots_dir
                )
                saved_jump_rate_plots.append(pr)

                # Save data CSV for local_jump
                data_jump_path = os.path.join(
                    data_dir,
                    f"{P}_{delta}_Cond1_step{step:04d}_jump_data.csv"
                )
                pd.DataFrame({
                    "arc_length": arclen,
                    "local_jump": local_jump
                }).to_csv(data_jump_path, index=False)
                saved_jump_data.append(data_jump_path)

                # Save data CSV for local_jump_rate
                data_rate_path = os.path.join(
                    data_dir,
                    f"{P}_{delta}_Cond1_step{step:04d}_jump_rate_data.csv"
                )
                pd.DataFrame({
                    "arc_length": arclen,
                    "local_jump_rate": local_jump_rate
                }).to_csv(data_rate_path, index=False)
                saved_rate_data.append(data_rate_path)

        # -----------------------------------
        # 6) Check Condition 2: x ≥ x_thresh
        # -----------------------------------
        if cond2_step is None:
            mask = df['x'] >= x_thresh
            if mask.any():
                region_rates = local_jump_rate[mask]
                max_rate2 = region_rates.max()
                if max_rate2 > rate_thresh:
                    idx2 = region_rates.idxmax()
                    x_at_idx2 = df['x'].iloc[idx2]
                    cond2_step = step
                    if verbose:
                        print(f"[{P}/{delta}] Cond2 met @ step {step:04d}, max_rate = {max_rate2:.2f} m/s at x = {x_at_idx2:.6f}")
                    # Plot local_jump
                    pj = plot_jump_and_save(
                        arclen=arclen,
                        local_jump_or_rate=local_jump,
                        highlight_idx=idx2,
                        step=step,
                        P=P,
                        delta=delta,
                        condition_label="Cond2",
                        flip_arc=flip_arc,
                        y_label="jump (m)",
                        curve_label="tangent jump",
                        plots_dir=plots_dir
                    )
                    saved_jump_plots.append(pj)
                    # Plot local_jump_rate
                    pr = plot_jump_and_save(
                        arclen=arclen,
                        local_jump_or_rate=local_jump_rate,
                        highlight_idx=idx2,
                        step=step,
                        P=P,
                        delta=delta,
                        condition_label="Cond2",
                        flip_arc=flip_arc,
                        y_label="jump_rate (m/s)",
                        curve_label="tangent jump rate",
                        plots_dir=plots_dir
                    )
                    saved_jump_rate_plots.append(pr)

                    # Save data CSV for local_jump
                    data_jump_path = os.path.join(
                        data_dir,
                        f"{P}_{delta}_Cond2_step{step:04d}_jump_data.csv"
                    )
                    pd.DataFrame({
                        "arc_length": arclen,
                        "local_jump": local_jump
                    }).to_csv(data_jump_path, index=False)
                    saved_jump_data.append(data_jump_path)

                    # Save data CSV for local_jump_rate
                    data_rate_path = os.path.join(
                        data_dir,
                        f"{P}_{delta}_Cond2_step{step:04d}_jump_rate_data.csv"
                    )
                    pd.DataFrame({
                        "arc_length": arclen,
                        "local_jump_rate": local_jump_rate
                    }).to_csv(data_rate_path, index=False)
                    saved_rate_data.append(data_rate_path)

        # If both conditions have been found, we can break early
        if cond1_step is not None and cond2_step is not None:
            break

    return (
        cond1_step, cond2_step,
        saved_jump_plots, saved_jump_rate_plots,
        saved_jump_data, saved_rate_data
    )


def main():
    args = parse_args()

    # Convert rotation angle from degrees to radians
    rot_rad = np.radians(args.rot_deg)

    # Root directory is current working directory (assumes P_cases are subfolders here)
    root_dir = "."

    # Ensure plot and data directories exist
    ensure_directory(args.plots_dir)
    data_dir = os.path.join(args.plots_dir, "data")
    ensure_directory(data_dir)

    summary_data = []
    all_jump_plots = []
    all_jump_rate_plots = []
    all_jump_data = []
    all_rate_data = []

    for P in args.P_cases:
        for delta in args.delta_cases:
            c1, c2, jp, jrp, jd, rd = scan_case(
                root_dir=root_dir,
                P=P,
                delta=delta,
                stem=args.stem,
                x0=args.x0,
                y0=args.y0,
                x_thresh=args.x_thresh,
                rate_thresh=args.rate_thresh,
                flip_index=args.flip_index,
                rot_rad=rot_rad,
                start_step=args.start_step,
                end_step=args.end_step,
                step_inc=args.step_increment,
                plots_dir=args.plots_dir,
                data_dir=data_dir,
                verbose=args.verbose
            )
            rel = None
            if c1 is not None and c2 is not None:
                rel = c2 - c1

            summary_data.append({
                "P_case": P,
                "delta_case": delta,
                "cond1_step": c1,
                "cond2_step": c2,
                "relative_step": rel
            })
            all_jump_plots.extend(jp)
            all_jump_rate_plots.extend(jrp)
            all_jump_data.extend(jd)
            all_rate_data.extend(rd)

    # Save summary CSV
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(args.output_csv, index=False)
    print(f"\nWrote summary to '{args.output_csv}'")

    if args.verbose:
        print("\nSaved jump plots:")
        for p in all_jump_plots:
            print("  " + p)
        print("\nSaved jump rate plots:")
        for p in all_jump_rate_plots:
            print("  " + p)
        print("\nSaved jump data CSVs:")
        for d in all_jump_data:
            print("  " + d)
        print("\nSaved jump rate data CSVs:")
        for d in all_rate_data:
            print("  " + d)


if __name__ == "__main__":
    main()