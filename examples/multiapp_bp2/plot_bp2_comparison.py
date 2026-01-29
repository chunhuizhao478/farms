#!/usr/bin/env python3
"""
Plot BP2 benchmark comparison: slip rate, slip, and shear stress vs time.
Compares MOOSE simulation results with SCEC SEAS benchmark data.
"""

import os

import matplotlib.pyplot as plt
import numpy as np

# File paths
script_dir = os.path.dirname(os.path.abspath(__file__))
moose_csv = os.path.join(script_dir, "main_elasticity_csv.csv")
benchmark_file = "./benchmark_data/bp2-qd-z0km-res.txt"


def load_moose_data(filename):
    """Load MOOSE CSV output."""
    data = np.genfromtxt(filename, delimiter=",", names=True)
    return {
        "time": data["time"],
        "slip": data["slip_z0"],
        "slip_rate": data["slip_rate_z0"],
        "shear_stress": data["shear_stress_z0"] / 1e6,  # Convert to MPa
        "state": data["state_z0"],
        "time_years": data["time_years"],
    }


def load_benchmark_data(filename):
    """Load SCEC SEAS benchmark data."""
    data = np.genfromtxt(
        filename,
        skip_header=16,
        names=["t", "slip", "slip_rate_log", "shear_stress", "state_log"],
    )
    return {
        "time": data["t"],
        "slip": data["slip"],
        "slip_rate": 10 ** data["slip_rate_log"],  # Convert from log10
        "shear_stress": data["shear_stress"],  # Already in MPa
        "state": 10 ** data["state_log"],  # Convert from log10
        "time_years": data["t"] / (365.25 * 24 * 3600),
    }


def plot_comparison(moose, benchmark, save_prefix="bp2_comparison"):
    """Create comparison plots."""

    # Convert time to years for plotting
    moose_years = moose["time"] / (365.25 * 24 * 3600)
    bench_years = benchmark["time"] / (365.25 * 24 * 3600)

    fig, axes = plt.subplots(4, 1, figsize=(12, 12), sharex=True)

    # --- Plot 1: Slip Rate (log scale) ---
    ax1 = axes[0]
    ax1.semilogy(
        moose_years, moose["slip_rate"], "b-", linewidth=1.5, label="MOOSE (this work)"
    )
    ax1.semilogy(
        bench_years,
        benchmark["slip_rate"],
        "r--",
        linewidth=1.5,
        label="Benchmark (Barbot)",
    )
    ax1.set_ylabel("Slip Rate (m/s)", fontsize=12)
    ax1.set_ylim([1e-15, 1e1])
    ax1.axhline(
        y=1e-3, color="gray", linestyle=":", alpha=0.5, label="Seismic threshold"
    )
    ax1.legend(loc="upper right", fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.set_title("BP2-QD Benchmark Comparison at z = 0 km (Free Surface)", fontsize=14)

    # --- Plot 2: Slip ---
    ax2 = axes[1]
    ax2.plot(moose_years, moose["slip"], "b-", linewidth=1.5, label="MOOSE")
    ax2.plot(bench_years, benchmark["slip"], "r--", linewidth=1.5, label="Benchmark")
    ax2.set_ylabel("Slip (m)", fontsize=12)
    ax2.legend(loc="upper left", fontsize=10)
    ax2.grid(True, alpha=0.3)

    # --- Plot 3: Shear Stress ---
    ax3 = axes[2]
    ax3.plot(moose_years, moose["shear_stress"], "b-", linewidth=1.5, label="MOOSE")
    ax3.plot(
        bench_years, benchmark["shear_stress"], "r--", linewidth=1.5, label="Benchmark"
    )
    ax3.set_ylabel("Shear Stress (MPa)", fontsize=12)
    ax3.ticklabel_format(
        axis="y", useOffset=False, style="plain"
    )  # Fix offset notation
    ax3.legend(loc="upper left", fontsize=10)
    ax3.grid(True, alpha=0.3)

    # --- Plot 4: State Variable (log scale) ---
    ax4 = axes[3]
    ax4.semilogy(moose_years, moose["state"], "b-", linewidth=1.5, label="MOOSE")
    ax4.semilogy(
        bench_years, benchmark["state"], "r--", linewidth=1.5, label="Benchmark"
    )
    ax4.set_ylabel("State Variable $\\theta$ (s)", fontsize=12)
    ax4.set_xlabel("Time (years)", fontsize=12)
    ax4.legend(loc="upper left", fontsize=10)
    ax4.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f"{save_prefix}_full.png", dpi=150, bbox_inches="tight")
    print(f"Saved: {save_prefix}_full.png")

    # --- Zoomed view of early time ---
    fig2, axes2 = plt.subplots(4, 1, figsize=(12, 12), sharex=True)

    # Zoom to first 0.1 years (about 36 days)
    zoom_limit = 0.1

    moose_mask = moose_years < zoom_limit
    bench_mask = bench_years < zoom_limit

    ax1 = axes2[0]
    ax1.semilogy(
        moose_years[moose_mask],
        moose["slip_rate"][moose_mask],
        "b-",
        linewidth=1.5,
        label="MOOSE",
    )
    ax1.semilogy(
        bench_years[bench_mask],
        benchmark["slip_rate"][bench_mask],
        "r--",
        linewidth=1.5,
        label="Benchmark",
    )
    ax1.set_ylabel("Slip Rate (m/s)", fontsize=12)
    ax1.legend(loc="upper right", fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.set_title(
        f"BP2-QD Early Time Comparison (first {zoom_limit} years)", fontsize=14
    )

    ax2 = axes2[1]
    ax2.plot(
        moose_years[moose_mask],
        moose["slip"][moose_mask] * 1000,
        "b-",
        linewidth=1.5,
        label="MOOSE",
    )
    ax2.plot(
        bench_years[bench_mask],
        benchmark["slip"][bench_mask] * 1000,
        "r--",
        linewidth=1.5,
        label="Benchmark",
    )
    ax2.set_ylabel("Slip (mm)", fontsize=12)
    ax2.legend(loc="upper left", fontsize=10)
    ax2.grid(True, alpha=0.3)

    ax3 = axes2[2]
    ax3.plot(
        moose_years[moose_mask],
        moose["shear_stress"][moose_mask],
        "b-",
        linewidth=1.5,
        label="MOOSE",
    )
    ax3.plot(
        bench_years[bench_mask],
        benchmark["shear_stress"][bench_mask],
        "r--",
        linewidth=1.5,
        label="Benchmark",
    )
    ax3.set_ylabel("Shear Stress (MPa)", fontsize=12)
    ax3.ticklabel_format(
        axis="y", useOffset=False, style="plain"
    )  # Fix offset notation
    ax3.legend(loc="upper left", fontsize=10)
    ax3.grid(True, alpha=0.3)

    ax4 = axes2[3]
    ax4.semilogy(
        moose_years[moose_mask],
        moose["state"][moose_mask],
        "b-",
        linewidth=1.5,
        label="MOOSE",
    )
    ax4.semilogy(
        bench_years[bench_mask],
        benchmark["state"][bench_mask],
        "r--",
        linewidth=1.5,
        label="Benchmark",
    )
    ax4.set_ylabel("State Variable $\\theta$ (s)", fontsize=12)
    ax4.set_xlabel("Time (years)", fontsize=12)
    ax4.legend(loc="upper left", fontsize=10)
    ax4.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f"{save_prefix}_zoom.png", dpi=150, bbox_inches="tight")
    print(f"Saved: {save_prefix}_zoom.png")

    plt.show()


def print_comparison_stats(moose, benchmark):
    """Print comparison statistics."""
    print("\n" + "=" * 60)
    print("BP2 BENCHMARK COMPARISON STATISTICS")
    print("=" * 60)

    # Initial values
    print("\nInitial Conditions (t=0):")
    print(
        f"  MOOSE:     τ = {moose['shear_stress'][0]:.4f} MPa, V = {moose['slip_rate'][0]:.2e} m/s, θ = {moose['state'][0]:.1f} s"
    )
    print(
        f"  Benchmark: τ = {benchmark['shear_stress'][0]:.4f} MPa, V = {benchmark['slip_rate'][0]:.2e} m/s, θ = {benchmark['state'][0]:.1f} s"
    )

    # Find first earthquake in each dataset (V > 1 mm/s)
    moose_eq_idx = np.where(moose["slip_rate"] > 1e-3)[0]
    bench_eq_idx = np.where(benchmark["slip_rate"] > 1e-3)[0]

    if len(moose_eq_idx) > 0:
        t_eq_moose = moose["time"][moose_eq_idx[0]]
        print(f"\nFirst Earthquake (V > 1 mm/s):")
        print(
            f"  MOOSE:     t = {t_eq_moose:.2e} s ({t_eq_moose / (365.25 * 24 * 3600):.4f} years)"
        )

    if len(bench_eq_idx) > 0:
        t_eq_bench = benchmark["time"][bench_eq_idx[0]]
        print(
            f"  Benchmark: t = {t_eq_bench:.2e} s ({t_eq_bench / (365.25 * 24 * 3600):.2f} years)"
        )

    # Stress evolution rate comparison
    print("\nStress Evolution (early time):")
    # Find stress at t ~ 1000 s
    moose_idx_1k = np.argmin(np.abs(moose["time"] - 1000))
    bench_idx_1k = np.argmin(np.abs(benchmark["time"] - 1000))

    if moose_idx_1k > 0:
        stress_rate_moose = (
            (moose["shear_stress"][moose_idx_1k] - moose["shear_stress"][0])
            / moose["time"][moose_idx_1k]
            * 1e6
        )  # Pa/s
        print(f"  MOOSE stress rate (0-1000s): {stress_rate_moose:.4f} Pa/s")

    if bench_idx_1k > 0:
        stress_rate_bench = (
            (benchmark["shear_stress"][bench_idx_1k] - benchmark["shear_stress"][0])
            / benchmark["time"][bench_idx_1k]
            * 1e6
        )
        print(f"  Benchmark stress rate (0-1000s): {stress_rate_bench:.4f} Pa/s")

    print("=" * 60)


if __name__ == "__main__":
    print("Loading data...")

    # Check if files exist
    if not os.path.exists(moose_csv):
        print(f"Error: MOOSE output not found: {moose_csv}")
        print("Run the simulation first.")
        exit(1)

    if not os.path.exists(benchmark_file):
        print(f"Error: Benchmark data not found: {benchmark_file}")
        exit(1)

    # Load data
    moose = load_moose_data(moose_csv)
    benchmark = load_benchmark_data(benchmark_file)

    print(f"MOOSE data: {len(moose['time'])} points, t_max = {moose['time'][-1]:.2e} s")
    print(
        f"Benchmark data: {len(benchmark['time'])} points, t_max = {benchmark['time'][-1]:.2e} s"
    )

    # Print statistics
    print_comparison_stats(moose, benchmark)

    # Create plots
    print("\nGenerating plots...")
    plot_comparison(moose, benchmark)
    print("Done!")
