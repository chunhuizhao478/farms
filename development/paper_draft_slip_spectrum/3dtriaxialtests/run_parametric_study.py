#!/usr/bin/env python3
"""Optimized confining-pressure parametric sweep for 3D triaxial test.

Performance improvements:
- Efficient CSV streaming (tail-like reading)
- Reduced file I/O overhead
- Optimized peak detection
- Better memory management
"""

from __future__ import annotations

import csv
import math
import os
import signal
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List

# ----------------------------------------------------------------------------


def status(section: str, message: str) -> None:
    """Consistent status banner for console output."""
    print(f"[{section}] {message}")

# User-adjustable settings
# ----------------------------------------------------------------------------
EXECUTABLE = Path("../../../farms-opt")
INPUT_FILE = Path("dynamic_solve_main.i")
OUTPUT_DIR = Path("parametric_results")

# Confining pressure targets, in MPa
# CONFINING_PRESSURES_MPA = [0, 50, 100, 200, 400, 600, 800]
CONFINING_PRESSURES_MPA = [25, 75, 150]

# Name of the CSV column containing the axial stress/force (sigma_1 proxy)
CSV_QUANTITY_TO_USE = "Fy"

# Cross-sectional area (m^2) used to convert the CSV force to stress
# (radius = 0.027 m).
FORCE_TO_STRESS_AREA: float | None = math.pi * (0.027 ** 2)

# Relative drop used to decide that the peak has been reached (e.g. 0.5 %)
PEAK_DROP_REL_TOL = 0 #it only has one peak value

# Require the drop condition to be met for this many consecutive samples
PEAK_DROP_MIN_SAMPLES = 3

# Seconds between checks for new CSV data while the solver runs
CSV_POLL_INTERVAL = 1.0  # Increased from 0.5s to reduce I/O overhead

# Only evaluate the peak/drop condition every N samples
PEAK_CHECK_INTERVAL = 1

# Enable debug output capture (may cause blocking on some systems)
CAPTURE_PROCESS_OUTPUT = False  # Set to True for debugging

# Conda environment for MOOSE execution
MOOSE_CONDA_ENV = "moose"  # Name of the conda environment with MOOSE dependencies

# Extra arguments for farms-opt (e.g., parallel launch, additional -P flags)
MPI_RUNNER = ["mpirun", "-np", "8"]
EXTRA_EXEC_ARGS: List[str] = []

# Buffer size for reading CSV files (lines at a time)
CSV_READ_BUFFER_LINES = 1000

# ----------------------------------------------------------------------------

@dataclass
class CaseResult:
    pressure_mpa: float
    pressure_pa: float
    sigma1: float
    i1: float
    sqrt_j2: float
    csv_file: Path

    def as_row(self) -> List[str]:
        return [
            f"{self.pressure_mpa:.1f}",
            f"{self.pressure_pa:.6e}",
            f"{self.sigma1:.6e}",
            f"{self.i1:.6e}",
            f"{self.sqrt_j2:.6e}",
            str(self.csv_file),
        ]


def terminate_process(proc: subprocess.Popen) -> None:
    if proc.poll() is None:
        try:
            os.killpg(proc.pid, signal.SIGTERM)
        except ProcessLookupError:
            return
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            proc.wait()


def parse_csv_line(line: str, header: List[str], column_name: str) -> tuple[float, float]:
    """Parse a single CSV line and extract time and target column value."""
    values = line.strip().split(',')
    if len(values) != len(header):
        raise ValueError(f"CSV line has {len(values)} values but header has {len(header)} columns")
    
    row_dict = dict(zip(header, values))
    time_val = float(row_dict.get("time", "nan"))
    raw_val = float(row_dict[column_name])
    return time_val, raw_val


def get_conda_run_command(env_name: str, cmd: List[str]) -> List[str]:
    """Wrap a command to run in a specific conda environment.
    
    Returns a command that activates the conda environment and runs the target command.
    """
    # Use 'conda run' which handles environment activation properly
    return ["conda", "run", "-n", env_name, "--no-capture-output"] + cmd


def run_command_until_peak(cmd: List[str], csv_file: Path, column_name: str) -> None:
    """Run the simulation, stopping early if the monitored stress drops after a peak.
    
    Optimized version with efficient CSV streaming.
    """
    csv_file.unlink(missing_ok=True)
    status("RUN", "Launching simulation")
    
    # Choose output handling based on configuration
    if CAPTURE_PROCESS_OUTPUT:
        proc = subprocess.Popen(
            cmd, 
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            bufsize=1,
            text=True,
            start_new_session=True,
        )
    else:
        # Silence solver output while allowing process-group termination
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
    
    terminated_early = False

    peak_stress = -math.inf
    drop_count = 0
    samples_processed = 0
    header: List[str] | None = None
    last_time = 0.0

    try:
        status("MONITOR", f"Waiting for CSV '{csv_file.name}'")
        
        # Wait for CSV file to be created
        max_wait = 60  # seconds
        wait_start = time.time()
        while not csv_file.exists():
            # Check if process died
            if proc.poll() is not None:
                if CAPTURE_PROCESS_OUTPUT:
                    stdout, stderr = proc.communicate()
                    
                    # Save output to debug file
                    debug_file = csv_file.parent / f"{csv_file.stem}_debug.log"
                    with debug_file.open('w') as f:
                        f.write("=== STDOUT ===\n")
                        f.write(stdout if stdout else "(empty)\n")
                        f.write("\n=== STDERR ===\n")
                        f.write(stderr if stderr else "(empty)\n")
                    
                    status("ERROR", f"Process terminated before CSV created (exit code {proc.returncode})")
                    status("ERROR", f"Debug output saved to {debug_file}")
                    if stderr:
                        stderr_lines = stderr.splitlines()
                        status("ERROR", f"Last stderr: {stderr_lines[-3:] if len(stderr_lines) >= 3 else stderr_lines}")
                else:
                    status("ERROR", f"Process terminated before CSV created (exit code {proc.returncode})")
                    status("ERROR", "Set CAPTURE_PROCESS_OUTPUT=True to save debug logs")
                
                raise RuntimeError(
                    f"Process terminated before CSV file created (exit code {proc.returncode}). "
                    f"The simulation may have failed during initialization."
                )
            
            if time.time() - wait_start > max_wait:
                # Process still running but no CSV - likely output configuration issue
                status("ERROR", f"CSV file not created after {max_wait}s. Process still running.")
                terminate_process(proc)
                raise TimeoutError(
                    f"CSV file {csv_file} not created after {max_wait}s. "
                    f"Check MOOSE Outputs configuration."
                )
            time.sleep(0.1)
        
        # Open file once and keep it open for streaming
        with csv_file.open('r') as fh:
            # Read header
            header_line = fh.readline()
            if not header_line:
                raise ValueError(f"Empty CSV file: {csv_file}")
            header = [col.strip() for col in header_line.strip().split(',')]
            
            if column_name not in header:
                raise KeyError(f"Column '{column_name}' not in CSV. Available: {header}")
            
            status("MONITOR", f"Monitoring column '{column_name}'")
            
            # Stream CSV data as it's written
            while True:
                line = fh.readline()
                
                if not line:
                    # No new data, check if process is still running
                    if proc.poll() is not None:
                        # Process finished, read any remaining data
                        for remaining_line in fh:
                            if not remaining_line.strip():
                                continue
                            try:
                                time_val, raw_val = parse_csv_line(remaining_line, header, column_name)
                                stress = to_stress(raw_val)
                                peak_stress = max(peak_stress, stress)
                            except (ValueError, KeyError):
                                continue
                        break
                    
                    # Check if we should terminate early
                    if terminated_early:
                        break
                    
                    time.sleep(CSV_POLL_INTERVAL)
                    continue
                
                # Parse the line
                if not line.strip():
                    continue
                
                try:
                    time_val, raw_val = parse_csv_line(line, header, column_name)
                except (ValueError, KeyError) as exc:
                    status("WARNING", f"Failed to parse CSV line: {exc}")
                    continue
                
                stress = to_stress(raw_val)
                samples_processed += 1
                last_time = time_val
                
                # Update peak or evaluate drop condition
                if stress > peak_stress:
                    peak_stress = stress
                    drop_count = 0
                    status(
                        "MONITOR",
                        f"t={time_val:.3e}s stress={stress:.3e}Pa → new peak (continue)",
                    )
                    continue

                if samples_processed % PEAK_CHECK_INTERVAL != 0:
                    status(
                        "MONITOR",
                        f"t={time_val:.3e}s stress={stress:.3e}Pa current_peak={peak_stress:.3e}Pa (continue)",
                    )
                    continue

                drop_condition = peak_stress != -math.inf and stress < peak_stress * (1 - PEAK_DROP_REL_TOL)
                if drop_condition:
                    drop_count += 1
                    if drop_count >= PEAK_DROP_MIN_SAMPLES:
                        status(
                            "MONITOR",
                            f"t={time_val:.3e}s stress={stress:.3e}Pa below peak for {drop_count} samples → stopping",
                        )
                        terminate_process(proc)
                        terminated_early = True
                        break
                    status(
                        "MONITOR",
                        f"t={time_val:.3e}s stress={stress:.3e}Pa below peak ({drop_count}/{PEAK_DROP_MIN_SAMPLES}) (continue)",
                    )
                else:
                    drop_count = 0
                    status(
                        "MONITOR",
                        f"t={time_val:.3e}s stress={stress:.3e}Pa current_peak={peak_stress:.3e}Pa (continue)",
                    )
                    
    except KeyboardInterrupt:
        status("MONITOR", "Interrupted by user")
        terminate_process(proc)
        raise
    except Exception:
        # On any error, try to terminate cleanly
        terminate_process(proc)
        raise
    finally:
        # Final cleanup
        if proc.poll() is None:
            terminate_process(proc)
    
    status("RUN", f"Simulation finished (processed {samples_processed} samples)")


def ensure_executable_available(exe: Path) -> None:
    status("SETUP", f"Checking executable at {exe}")
    if not exe.exists():
        raise FileNotFoundError(
            f"Cannot locate executable '{exe}'. Adjust EXECUTABLE in the script."
        )
    if not exe.is_file():
        raise FileNotFoundError(f"'{exe}' is not a file. Check your build output.")
    status("SETUP", "Executable available")


def check_conda_environment(env_name: str) -> None:
    """Verify that the specified conda environment exists."""
    status("SETUP", f"Checking conda environment '{env_name}'")
    try:
        result = subprocess.run(
            ["conda", "env", "list"],
            capture_output=True,
            text=True,
            timeout=10
        )
        if result.returncode != 0:
            raise RuntimeError("Failed to list conda environments")
        
        # Check if environment name appears in the output
        if env_name not in result.stdout:
            raise EnvironmentError(
                f"Conda environment '{env_name}' not found. "
                f"Please create it or update MOOSE_CONDA_ENV in the script."
            )
        status("SETUP", f"Conda environment '{env_name}' found")
    except FileNotFoundError:
        raise RuntimeError(
            "conda command not found. Ensure conda is installed and in PATH."
        )


def load_peak_value(csv_path: Path, column_name: str) -> float:
    """Return the maximum value in the specified CSV column.
    
    Optimized to use csv module efficiently.
    """
    peak = -math.inf
    with csv_path.open('r') as fh:
        reader = csv.DictReader(fh)
        if column_name not in reader.fieldnames:
            raise KeyError(
                f"Column '{column_name}' not found in {csv_path}. "
                f"Available columns: {reader.fieldnames}"
            )
        
        for row in reader:
            try:
                value = float(row[column_name])
                peak = max(peak, value)
            except (ValueError, KeyError) as exc:
                status("WARNING", f"Skipping invalid row: {exc}")
                continue
        
        if peak == -math.inf:
            raise ValueError(f"No valid data rows found in {csv_path}")
    
    return peak


def to_stress(value: float) -> float:
    if FORCE_TO_STRESS_AREA is None:
        return value
    if FORCE_TO_STRESS_AREA <= 0:
        raise ValueError("FORCE_TO_STRESS_AREA must be positive if provided")
    return value / FORCE_TO_STRESS_AREA


def compute_invariants(sigma1: float, sigma_conf: float) -> tuple[float, float]:
    i1 = sigma1 + 2.0 * sigma_conf
    sqrt_j2 = abs(sigma1 - sigma_conf) / math.sqrt(3.0)
    return i1, sqrt_j2


def build_case_name(pressure_mpa: float) -> str:
    if pressure_mpa.is_integer():
        return f"sigma_{int(pressure_mpa):03d}MPa"
    return f"sigma_{pressure_mpa:06.2f}MPa"


def run_case(pressure_mpa: float) -> CaseResult:
    pressure_pa = pressure_mpa * 1e6
    case_name = build_case_name(pressure_mpa)
    csv_base = f"{case_name}_main"
    csv_base_path = OUTPUT_DIR / csv_base
    csv_file = csv_base_path.with_suffix(".csv")

    # Build the base command (without conda wrapper yet)
    base_cmd: List[str] = MPI_RUNNER + [str(EXECUTABLE), "-i", str(INPUT_FILE)]
    base_cmd += ["-P", f"confining_pressure={pressure_pa}"]

    # Ensure each output object writes to a unique file base to avoid overwrites
    base_cmd += ["-P", f"Outputs/csv/file_base={csv_base_path}"]
    base_cmd += ["-P", f"Outputs/exodus/file_base={OUTPUT_DIR / case_name}"]
    base_cmd += EXTRA_EXEC_ARGS
    
    # Wrap with conda environment activation
    cmd = get_conda_run_command(MOOSE_CONDA_ENV, base_cmd)

    status("CASE", f"=== Running case {case_name} (sigma_c = {pressure_mpa} MPa) ===")
    status("CASE", f"Command: {' '.join(base_cmd)}")
    status("CASE", f"Conda environment: {MOOSE_CONDA_ENV}")
    
    case_start = time.time()
    run_command_until_peak(cmd, csv_file, CSV_QUANTITY_TO_USE)
    case_duration = time.time() - case_start

    if not csv_file.exists():
        raise FileNotFoundError(
            f"Expected CSV '{csv_file}' was not created. Check CSV output settings."
        )

    peak_value = load_peak_value(csv_file, CSV_QUANTITY_TO_USE)
    sigma1 = to_stress(peak_value)
    i1, sqrt_j2 = compute_invariants(sigma1, pressure_pa)

    status(
        "CASE",
        f"Completed {case_name} in {case_duration:.1f}s: "
        f"sigma1_peak = {sigma1:.3e} Pa, I1 = {i1:.3e} Pa, sqrt(J2) = {sqrt_j2:.3e} Pa"
    )

    return CaseResult(
        pressure_mpa=pressure_mpa,
        pressure_pa=pressure_pa,
        sigma1=sigma1,
        i1=i1,
        sqrt_j2=sqrt_j2,
        csv_file=csv_file,
    )


def write_summary(results: Iterable[CaseResult]) -> Path:
    summary_path = OUTPUT_DIR / "summary.csv"
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    with summary_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            "sigma_c_MPa", "sigma_c_Pa", "sigma1_peak_Pa", 
            "I1_Pa", "sqrtJ2_Pa", "csv_file"
        ])
        for result in results:
            writer.writerow(result.as_row())

    status("SUMMARY", f"Summary written to {summary_path}")
    return summary_path


def main() -> None:
    start_time = time.time()
    status("SETUP", "Initializing parametric study")
    
    # Check conda environment first
    check_conda_environment(MOOSE_CONDA_ENV)
    
    ensure_executable_available(EXECUTABLE)
    if not INPUT_FILE.exists():
        raise FileNotFoundError(f"Input file '{INPUT_FILE}' not found")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    status("SETUP", f"Output directory: {OUTPUT_DIR}")

    for path in OUTPUT_DIR.glob("*.csv"):
        path.unlink(missing_ok=True)
    for path in OUTPUT_DIR.glob("*.e"):
        path.unlink(missing_ok=True)
    for path in OUTPUT_DIR.glob("*.e.*"):
        path.unlink(missing_ok=True)
    status("SETUP", "Cleared previous CSV/Exodus files")

    summary_path = OUTPUT_DIR / "summary.csv"
    with summary_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            "sigma_c_MPa",
            "sigma_c_Pa",
            "sigma1_peak_Pa",
            "I1_Pa",
            "sqrtJ2_Pa",
            "csv_file",
        ])
    status("SETUP", f"Initialized summary file at {summary_path}")

    results: List[CaseResult] = []
    total_cases = len(CONFINING_PRESSURES_MPA)
    for idx, pressure_mpa in enumerate(CONFINING_PRESSURES_MPA, start=1):
        status("PROGRESS", f"Running case {idx}/{total_cases}")
        result = run_case(float(pressure_mpa))
        results.append(result)
        with summary_path.open("a", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(result.as_row())
        status("PROGRESS", f"Recorded result for {result.pressure_mpa:.1f} MPa")

    summary_path = write_summary(results)

    i1_list = [r.i1 for r in results]
    sqrt_j2_list = [r.sqrt_j2 for r in results]
    status("SUMMARY", f"I1 values (Pa): {i1_list}")
    status("SUMMARY", f"sqrt(J2) values (Pa): {sqrt_j2_list}")
    status("SUMMARY", f"Results saved to {summary_path}")
    
    total_time = time.time() - start_time
    status("SUMMARY", f"Total execution time: {total_time:.1f}s ({total_time/60:.1f} min)")


if __name__ == "__main__":
    main()
