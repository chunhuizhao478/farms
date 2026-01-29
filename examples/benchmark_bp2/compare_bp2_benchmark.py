#!/usr/bin/env python3
"""
Compare MOOSE SEAS BP2 simulation results with SCEC benchmark data.

This script reads:
1. Benchmark data from SCEC (bp2-qd-z0km-res.txt)
2. MOOSE simulation output (bp2_800m_z0km.csv)

And produces comparison plots of:
- Slip vs time
- Slip rate (log10) vs time
- Shear stress vs time
- State variable (log10) vs time

Usage:
    python compare_bp2_benchmark.py [benchmark_file] [simulation_file]
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os

# Seconds per year
SECONDS_PER_YEAR = 365.25 * 24 * 3600


def read_benchmark_data(filename):
    """
    Read SCEC benchmark data file.

    Format:
    - Header lines starting with #
    - Column names: t slip slip_rate shear_stress state
    - Data columns: time(s), slip(m), log10(V), stress(MPa), log10(theta)
    """
    # Skip header lines and read data
    data = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.strip().startswith('t '):
                continue  # Skip column header line
            parts = line.split()
            if len(parts) >= 5:
                try:
                    row = [float(x) for x in parts[:5]]
                    data.append(row)
                except ValueError:
                    continue

    data = np.array(data)

    return {
        'time': data[:, 0],           # seconds
        'slip': data[:, 1],           # meters
        'slip_rate_log': data[:, 2],  # log10(m/s)
        'shear_stress': data[:, 3],   # MPa
        'state_log': data[:, 4]       # log10(s)
    }


def read_moose_csv(filename):
    """
    Read MOOSE CSV output file.

    Expected columns:
    - time
    - slip_z0
    - slip_rate_z0
    - shear_stress_z0
    - state_z0
    """
    df = pd.read_csv(filename)

    # Handle different column naming conventions
    result = {'time': df['time'].values}

    if 'slip_z0' in df.columns:
        result['slip'] = df['slip_z0'].values
    elif 'slip' in df.columns:
        result['slip'] = df['slip'].values

    if 'slip_rate_z0' in df.columns:
        result['slip_rate'] = df['slip_rate_z0'].values
    elif 'slip_rate' in df.columns:
        result['slip_rate'] = df['slip_rate'].values

    if 'shear_stress_z0' in df.columns:
        result['shear_stress'] = df['shear_stress_z0'].values / 1e6  # Pa to MPa
    elif 'traction' in df.columns:
        result['shear_stress'] = df['traction'].values / 1e6

    if 'state_z0' in df.columns:
        result['state'] = df['state_z0'].values
    elif 'state_variable' in df.columns:
        result['state'] = df['state_variable'].values

    # Compute log10 values
    if 'slip_rate' in result:
        result['slip_rate_log'] = np.log10(np.maximum(result['slip_rate'], 1e-20))
    if 'state' in result:
        result['state_log'] = np.log10(np.maximum(result['state'], 1e-20))

    return result


def plot_comparison(benchmark, simulation, output_prefix='bp2_comparison'):
    """
    Create comparison plots.
    """
    # Convert times to years for plotting
    bench_time_yr = benchmark['time'] / SECONDS_PER_YEAR
    sim_time_yr = simulation['time'] / SECONDS_PER_YEAR

    # Create figure with 4 subplots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # 1. Slip vs time
    ax = axes[0, 0]
    ax.plot(bench_time_yr, benchmark['slip'], 'b-', label='Benchmark', linewidth=1)
    if 'slip' in simulation:
        ax.plot(sim_time_yr, simulation['slip'], 'r--', label='MOOSE', linewidth=1)
    ax.set_xlabel('Time (years)')
    ax.set_ylabel('Slip (m)')
    ax.set_title('Slip at z = 0 km')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 2. Slip rate (log10) vs time
    ax = axes[0, 1]
    ax.plot(bench_time_yr, benchmark['slip_rate_log'], 'b-', label='Benchmark', linewidth=1)
    if 'slip_rate_log' in simulation:
        ax.plot(sim_time_yr, simulation['slip_rate_log'], 'r--', label='MOOSE', linewidth=1)
    ax.set_xlabel('Time (years)')
    ax.set_ylabel('log₁₀(Slip rate) (log₁₀ m/s)')
    ax.set_title('Slip Rate at z = 0 km')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.axhline(y=-3, color='k', linestyle=':', alpha=0.5, label='V = 1 mm/s (seismic)')

    # 3. Shear stress vs time
    ax = axes[1, 0]
    ax.plot(bench_time_yr, benchmark['shear_stress'], 'b-', label='Benchmark', linewidth=1)
    if 'shear_stress' in simulation:
        ax.plot(sim_time_yr, simulation['shear_stress'], 'r--', label='MOOSE', linewidth=1)
    ax.set_xlabel('Time (years)')
    ax.set_ylabel('Shear Stress (MPa)')
    ax.set_title('Shear Stress at z = 0 km')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 4. State variable (log10) vs time
    ax = axes[1, 1]
    ax.plot(bench_time_yr, benchmark['state_log'], 'b-', label='Benchmark', linewidth=1)
    if 'state_log' in simulation:
        ax.plot(sim_time_yr, simulation['state_log'], 'r--', label='MOOSE', linewidth=1)
    ax.set_xlabel('Time (years)')
    ax.set_ylabel('log₁₀(State) (log₁₀ s)')
    ax.set_title('State Variable at z = 0 km')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'{output_prefix}.png', dpi=150)
    plt.savefig(f'{output_prefix}.pdf')
    print(f"Saved comparison plots to {output_prefix}.png and {output_prefix}.pdf")

    return fig


def compute_errors(benchmark, simulation, time_range=None):
    """
    Compute error metrics between benchmark and simulation.

    Args:
        benchmark: Dictionary with benchmark data
        simulation: Dictionary with simulation data
        time_range: Optional tuple (t_min, t_max) in seconds

    Returns:
        Dictionary with error metrics
    """
    # Interpolate simulation data to benchmark times
    from scipy import interpolate

    bench_time = benchmark['time']

    if time_range:
        mask = (bench_time >= time_range[0]) & (bench_time <= time_range[1])
        bench_time = bench_time[mask]

    errors = {}

    for var in ['slip', 'slip_rate_log', 'shear_stress', 'state_log']:
        if var not in simulation:
            continue

        # Interpolate simulation to benchmark times
        sim_interp = interpolate.interp1d(
            simulation['time'],
            simulation[var if var != 'slip_rate_log' else 'slip_rate_log'],
            bounds_error=False,
            fill_value='extrapolate'
        )

        sim_at_bench = sim_interp(bench_time)

        if time_range:
            bench_var = benchmark[var][mask]
        else:
            bench_var = benchmark[var]

        # L2 relative error
        l2_error = np.sqrt(np.mean((sim_at_bench - bench_var)**2))
        l2_rel = l2_error / np.sqrt(np.mean(bench_var**2)) if np.any(bench_var != 0) else np.inf

        # Maximum absolute error
        max_error = np.max(np.abs(sim_at_bench - bench_var))

        errors[var] = {
            'l2_error': l2_error,
            'l2_relative': l2_rel,
            'max_error': max_error
        }

    return errors


def print_summary(benchmark, simulation, errors=None):
    """
    Print summary statistics.
    """
    print("\n" + "="*60)
    print("SEAS BP2 Benchmark Comparison Summary")
    print("="*60)

    print(f"\nBenchmark data: {len(benchmark['time'])} time steps")
    print(f"  Time range: {benchmark['time'][0]:.2e} - {benchmark['time'][-1]:.2e} s")
    print(f"            = {benchmark['time'][0]/SECONDS_PER_YEAR:.2f} - {benchmark['time'][-1]/SECONDS_PER_YEAR:.2f} years")

    print(f"\nSimulation data: {len(simulation['time'])} time steps")
    print(f"  Time range: {simulation['time'][0]:.2e} - {simulation['time'][-1]:.2e} s")
    print(f"            = {simulation['time'][0]/SECONDS_PER_YEAR:.2f} - {simulation['time'][-1]/SECONDS_PER_YEAR:.2f} years")

    if errors:
        print("\nError Metrics (over common time range):")
        print("-"*60)
        print(f"{'Variable':<20} {'L2 Error':<15} {'L2 Relative':<15} {'Max Error':<15}")
        print("-"*60)
        for var, err in errors.items():
            print(f"{var:<20} {err['l2_error']:<15.4e} {err['l2_relative']:<15.4e} {err['max_error']:<15.4e}")

    print("="*60)


def main():
    # Default file paths
    benchmark_file = "/Users/chunhuizhao/Desktop/bp2-qd-z0km-res.txt"
    simulation_file = "bp2_800m_z0km.csv"

    # Override with command line arguments
    if len(sys.argv) >= 2:
        benchmark_file = sys.argv[1]
    if len(sys.argv) >= 3:
        simulation_file = sys.argv[2]

    print(f"Reading benchmark data from: {benchmark_file}")
    print(f"Reading simulation data from: {simulation_file}")

    # Read data
    try:
        benchmark = read_benchmark_data(benchmark_file)
        print(f"  Loaded {len(benchmark['time'])} benchmark data points")
    except Exception as e:
        print(f"Error reading benchmark file: {e}")
        return 1

    try:
        simulation = read_moose_csv(simulation_file)
        print(f"  Loaded {len(simulation['time'])} simulation data points")
    except FileNotFoundError:
        print(f"Simulation file not found: {simulation_file}")
        print("Creating comparison with benchmark data only...")
        simulation = {'time': np.array([0])}
    except Exception as e:
        print(f"Error reading simulation file: {e}")
        simulation = {'time': np.array([0])}

    # Compute errors if both datasets available
    errors = None
    if len(simulation['time']) > 1:
        # Find common time range
        t_min = max(benchmark['time'][0], simulation['time'][0])
        t_max = min(benchmark['time'][-1], simulation['time'][-1])

        if t_max > t_min:
            try:
                errors = compute_errors(benchmark, simulation, (t_min, t_max))
            except Exception as e:
                print(f"Warning: Could not compute errors: {e}")

    # Print summary
    print_summary(benchmark, simulation, errors)

    # Create plots
    output_prefix = os.path.splitext(simulation_file)[0] + '_comparison'
    try:
        plot_comparison(benchmark, simulation, output_prefix)
    except Exception as e:
        print(f"Error creating plots: {e}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
