#!/usr/bin/env python3
"""Fit linear Drucker–Prager envelopes for multiple summary files.

For each summary CSV, perform sqrt(J2) = A + B I1 regression, plot the data
with the fitted line, and compute cohesion c and friction angle phi using two
alternative parameterizations.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Iterable, Tuple

import matplotlib.pyplot as plt
import numpy as np

SUMMARY_FILES = [
    (Path("summary_eps1em7.csv"), "$\hat{\dot{\epsilon}} = 10^{-7}$ 1/s, $\hat{C_d} = 10$ 1/s"),
    (Path("summary_eps1em8.csv"), "$\hat{\dot{\epsilon}} = 10^{-8}$ 1/s, $\hat{C_d} = 10$ 1/s"),
    (Path("summary_eps1em9.csv"), "$\hat{\dot{\epsilon}} = 10^{-9}$ 1/s, $\hat{C_d} = 10$ 1/s"),
]
DEFAULT_COLORS = ["C0", "C1", "C2", "C3"]
OUTPUT_FIG = Path("drucker_prager_fit.png")


def load_summary(path: Path) -> Tuple[np.ndarray, np.ndarray]:
    """Return arrays of I1 and sqrt(J2) from the summary CSV."""
    i1_vals = []
    j2_vals = []
    with path.open("r", newline="") as fh:
        reader = csv.DictReader(fh)
        if "I1_Pa" not in reader.fieldnames or "sqrtJ2_Pa" not in reader.fieldnames:
            raise KeyError("summary CSV must contain 'I1_Pa' and 'sqrtJ2_Pa' columns")
        for row in reader:
            try:
                i1_vals.append(float(row["I1_Pa"]))
                j2_vals.append(float(row["sqrtJ2_Pa"]))
            except ValueError as exc:
                raise ValueError(f"Invalid numeric entry in row {row}") from exc
    if not i1_vals:
        raise ValueError(f"No data rows found in summary file {path}")
    return np.array(i1_vals), np.array(j2_vals)


def fit_linear(i1: np.ndarray, sqrt_j2: np.ndarray) -> Tuple[float, float]:
    """Least-squares fit of sqrt(J2) = A + B * I1."""
    coeffs = np.polyfit(i1, sqrt_j2, 1)
    B, A = coeffs[0], coeffs[1]
    return float(A), float(B)


def solve_cohesion_friction(A: float, B: float) -> Tuple[float, float]:
    """Solve for cohesion c and friction angle phi (radians) using set 1."""
    sqrt3 = math.sqrt(3.0)
    if abs(B) < 1e-14:
        sin_phi = 0.0
    else:
        sin_phi = (3.0 * B * sqrt3) / (2.0 + B * sqrt3)
    sin_phi = max(-1.0, min(1.0, sin_phi))
    phi = math.asin(sin_phi)
    cos_phi = math.cos(phi)
    denom = sqrt3 * (3.0 - sin_phi)
    if abs(denom) < 1e-14:
        raise ZeroDivisionError("Denominator nearly zero while solving for cohesion")
    c = (A * denom) / (6.0 * cos_phi)
    return c, phi


def solve_cohesion_friction_alt(A: float, B: float) -> Tuple[float, float]:
    """Solve for cohesion c and phi using alternative parameterization."""
    denominator = 1.0 - 3.0 * B * B
    if abs(denominator) < 1e-14:
        raise ZeroDivisionError("Cannot solve for sin(phi): denominator near zero")
    sin_phi_sq = 9.0 * B * B / denominator
    sin_phi_sq = max(0.0, sin_phi_sq)
    sin_phi = math.sqrt(sin_phi_sq)
    sin_phi = max(-1.0, min(1.0, sin_phi))
    phi = math.asin(sin_phi)
    cos_phi = math.cos(phi)
    denom = math.sqrt(9.0 + 3.0 * sin_phi * sin_phi)
    if abs(denom) < 1e-14:
        raise ZeroDivisionError("Cannot compute cohesion: denominator near zero")
    c = (A * denom) / (3.0 * cos_phi)
    return c, phi


def main() -> None:
    plt.figure(figsize=(6, 4))
    summaries: list[Tuple[str, float, float, float, float, float, float]] = []

    for idx, (path, user_label) in enumerate(SUMMARY_FILES):
        if not path.exists():
            print(f"Warning: {path} not found; skipping")
            continue
        color = DEFAULT_COLORS[idx % len(DEFAULT_COLORS)]
        label = user_label or path.stem

        i1, sqrt_j2 = load_summary(path)
        A, B = fit_linear(i1, sqrt_j2)
        c1, phi1 = solve_cohesion_friction(A, B)
        c2, phi2 = solve_cohesion_friction_alt(A, B)
        summaries.append((label, A, B, c1, phi1, c2, phi2))

        i1_grid = np.linspace(i1.min(), i1.max(), 200)
        sqrt_j2_fit = A + B * i1_grid

        plt.scatter(i1/1e6, sqrt_j2/1e6, color=color, alpha=0.6, label=f"{label} data")
        plt.plot(i1_grid/1e6, sqrt_j2_fit/1e6, color=color, linestyle="--", label=f"{label} fit")

    plt.xlabel(r"$\sigma_m$ (MPa)")
    plt.ylabel(r"$\sqrt{J_2}$ (MPa)")
    plt.title("Drucker–Prager Yield Criterion for Different $\hat{\dot{\epsilon}}$")
    plt.legend(fontsize=7)
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    plt.savefig(OUTPUT_FIG, dpi=200)
    print(f"Saved plot to {OUTPUT_FIG}")

    for label, A, B, c1, phi1, c2, phi2 in summaries:
        print(f"\nSummary for {label}:")
        print(f"  sqrt(J2) = {A:.3e} + ({B:.3e}) * I1")
        print("  --- Parameter Set 1 ---")
        print(f"    Cohesion c = {c1:.3e} Pa")
        print(f"    Friction angle phi = {math.degrees(phi1):.2f} degrees")
        print("  --- Parameter Set 2 ---")
        print(f"    Cohesion c = {c2:.3e} Pa")
        print(f"    Friction angle phi = {math.degrees(phi2):.2f} degrees")


if __name__ == "__main__":
    main()
