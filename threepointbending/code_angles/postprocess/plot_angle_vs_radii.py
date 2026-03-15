"""
Postprocessing script: Angle of Intersection with Hole vs Number of Radii

Reproduces the experimental figure and overlays simulation results.

Angle convention (looking at the circle from outside):
  0°   = bottom of hole
  90°  = right side
  -90° = left side
  ±180° = top of hole

Number of Radii = dx / hole_radius (x-offset from notch tip divided by hole radius)

Usage:
  python plot_angle_vs_radii.py
"""

import os
import re
from math import atan2, degrees

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
# Paths
# =============================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
GEO_DIR = os.path.join(SCRIPT_DIR, "..", "..", "meshfile_angles")
OUTPUT_DIR = SCRIPT_DIR

# =============================================================================
# Hole parameters (constant across all cases)
# =============================================================================
NOTCH_TIP_X = 0.0140  # m
NOTCH_TIP_Y = 0.0016  # m
DY = 0.0016  # m (vertical offset, constant)
HOLE_RADIUS = 0.001  # m (1 mm)

# =============================================================================
# Experimental data (extracted from figure)
# Format: (number_of_radii, angle_degrees)
# =============================================================================

# 2mm Trial 1 (blue squares with dashed line)
exp_2mm_trial1 = np.array(
    [
        [-4.0, -180],
        [-3.0, -180],
        [-2.0, -70],
        [-1.0, -28],
        [-0.5, -40],
        [0.0, 15],
        [0.5, -35],
        [1.0, 20],
        [1.0, -40],
        [2.0, 65],
        [2.5, 90],
        [3.0, 175],
        [4.0, 180],
    ]
)

# 2mm Trial 2 (red squares)
exp_2mm_trial2 = np.array(
    [
        [-1.0, -45],
        [-0.5, -40],
        [0.0, -40],
        [0.5, -30],
        [1.0, 20],
    ]
)

# 1.5mm (green squares)
exp_1_5mm = np.array(
    [
        [0.0, 10],
        [0.5, 55],
        [1.0, 65],
        [1.5, 65],
    ]
)

# Dashed trend line (digitized from experimental figure curve.txt)
_curve_data = np.loadtxt(os.path.join(SCRIPT_DIR, "curve.txt"), delimiter=",")
trend_x = _curve_data[:, 0]
trend_y = _curve_data[:, 1]


# =============================================================================
# Geo file parser
# =============================================================================
def parse_geo_file(filepath):
    """Extract hole parameters from a .geo file."""
    params = {}
    with open(filepath) as f:
        for line in f:
            for key in ["dx", "dy", "hole_radius", "notch_tip_x", "notch_tip_y"]:
                m = re.match(rf"^{key}\s*=\s*([0-9.eE+-]+)\s*;", line)
                if m:
                    params[key] = float(m.group(1))
    params["hole_center_x"] = params["notch_tip_x"] + params["dx"]
    params["hole_center_y"] = params["notch_tip_y"] + params["dy"]
    return params


# =============================================================================
# Angle calculation
# =============================================================================
def compute_intersection_angle(hole_center, contact_point):
    """
    Compute angle of intersection using the convention:
      0°   = bottom of circle
      90°  = right
      -90° = left
      ±180° = top

    Parameters
    ----------
    hole_center : (cx, cy) in meters
    contact_point : (px, py) in meters

    Returns
    -------
    angle in degrees
    """
    cx, cy = hole_center
    px, py = contact_point
    dx = px - cx
    dy = py - cy
    angle = degrees(atan2(dx, -dy))
    return angle


# =============================================================================
# Build simulation cases from .geo files
# =============================================================================
def build_simulation_cases():
    """Read all .geo files and build case list."""
    labels = [
        "0.00",
        "0.25",
        "0.50",
        "0.75",
        "1.00",
        "1.25",
        "1.50",
        "1.75",
        "2.00",
        "2.25",
        "2.50",
        "2.75",
        "3.00",
    ]
    cases = []
    for label in labels:
        geo_file = os.path.join(GEO_DIR, f"mesh_whole_3d_x{label}mm_y1.6mm.geo")
        if not os.path.exists(geo_file):
            print(f"Warning: {geo_file} not found, skipping")
            continue
        params = parse_geo_file(geo_file)
        num_radii = params["dx"] / params["hole_radius"]
        cases.append(
            {
                "label": label,
                "dx": params["dx"],
                "num_radii": num_radii,
                "hole_center": (params["hole_center_x"], params["hole_center_y"]),
                "hole_radius": params["hole_radius"],
            }
        )
    return cases


# =============================================================================
# Simulation contact points (USER INPUT)
# =============================================================================
# For each case label, provide the (x, y) coordinate where the damage band
# contacts the hole circle, read from ParaView at z = 0.008 plane.
# Set to None for cases not yet measured.
#
# Example: 'x0.00mm': (0.0140, 0.0022)
#
CONTACT_POINTS = {
    "0.00": (0.0148089),
    "0.25": None,
    "0.50": None,
    "0.75": None,
    "1.00": None,
    "1.25": None,
    "1.50": None,
    "1.75": None,
    "2.00": None,
    "2.25": None,
    "2.50": None,
    "2.75": None,
    "3.00": None,
}


# =============================================================================
# Compute simulation angles
# =============================================================================
def compute_simulation_data(cases, contact_points):
    """Compute angles for all cases with contact points."""
    num_radii_list = []
    angle_list = []
    for case in cases:
        cp = contact_points.get(case["label"])
        if cp is not None:
            angle = compute_intersection_angle(case["hole_center"], cp)
            num_radii_list.append(case["num_radii"])
            angle_list.append(angle)
    return np.array(num_radii_list), np.array(angle_list)


# =============================================================================
# Draw angle convention inset
# =============================================================================
def draw_angle_convention(ax):
    """Draw the angle convention circle diagram as an inset."""
    inset = ax.inset_axes([0.02, 0.55, 0.22, 0.40])
    theta = np.linspace(0, 2 * np.pi, 100)
    inset.plot(np.cos(theta), np.sin(theta), "b-", linewidth=1.2)

    # Dashed crosshair
    inset.plot([0, 0], [-1.1, 1.1], "b--", linewidth=0.5, alpha=0.5)
    inset.plot([-1.1, 1.1], [0, 0], "b--", linewidth=0.5, alpha=0.5)

    # Arrows and labels
    inset.annotate(
        "",
        xy=(0.15, 0.85),
        xytext=(0, 0),
        arrowprops=dict(arrowstyle="->", color="b", lw=1.2),
    )
    fs = 8
    inset.text(0, -1.4, "0°", ha="center", va="top", fontsize=fs, color="b")
    inset.text(1.3, 0, "90°", ha="left", va="center", fontsize=fs, color="b")
    inset.text(-1.3, 0, "-90°", ha="right", va="center", fontsize=fs, color="b")
    inset.text(0, 1.35, "±180°", ha="center", va="bottom", fontsize=fs, color="b")
    inset.text(0.45, 0.45, "R", ha="center", va="center", fontsize=fs, color="b")

    inset.set_xlim(-1.8, 1.8)
    inset.set_ylim(-1.8, 1.8)
    inset.set_aspect("equal")
    inset.axis("off")


# =============================================================================
# Plotting
# =============================================================================
def create_plot(sim_num_radii, sim_angles, output_path):
    """Create the angle vs number of radii plot."""
    fig, ax = plt.subplots(figsize=(10, 5))

    # Trend line (experimental best fit)
    ax.plot(
        trend_x, trend_y, "k--", linewidth=1, alpha=0.5, label="Experiment (best fit)"
    )

    # Simulation data
    if len(sim_num_radii) > 0:
        ax.plot(
            sim_num_radii,
            sim_angles,
            "o",
            color="black",
            markersize=8,
            markerfacecolor="black",
            markeredgewidth=1.5,
            label="Simulation (1mm radius)",
        )

    # Axis settings
    ax.set_xlim(-4.5, 4.5)
    ax.set_ylim(-180, 180)
    ax.set_yticks([-180, -120, -60, 0, 60, 120, 180])
    ax.set_xticks(range(-4, 5))
    ax.set_xlabel("Number of Radii", fontsize=13)
    ax.set_ylabel("Angle of Intersection with Hole (degrees)", fontsize=13)

    # Grid
    ax.axhline(y=0, color="gray", linewidth=0.5, linestyle="-")
    ax.axvline(x=0, color="gray", linewidth=0.5, linestyle="-")
    ax.grid(True, axis="x", linestyle=":", alpha=0.4)
    ax.grid(True, axis="y", linestyle=":", alpha=0.3)

    # Legend
    ax.legend(loc="lower right", fontsize=10, framealpha=0.9)

    # Angle convention inset
    draw_angle_convention(ax)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.savefig(output_path.replace(".png", ".pdf"), bbox_inches="tight")
    print(f"Saved: {output_path}")
    print(f"Saved: {output_path.replace('.png', '.pdf')}")


# =============================================================================
# Main
# =============================================================================
if __name__ == "__main__":
    # Build cases from .geo files
    cases = build_simulation_cases()
    print(f"Found {len(cases)} simulation cases:")
    for c in cases:
        print(
            f"  x={c['label']}mm  num_radii={c['num_radii']:.2f}  "
            f"hole_center=({c['hole_center'][0]:.5f}, {c['hole_center'][1]:.5f})"
        )

    # Compute simulation angles
    sim_nr, sim_ang = compute_simulation_data(cases, CONTACT_POINTS)
    if len(sim_nr) > 0:
        print(f"\nSimulation results ({len(sim_nr)} points):")
        for nr, ang in zip(sim_nr, sim_ang):
            print(f"  num_radii={nr:.2f}  angle={ang:.1f}°")
    else:
        print("\nNo simulation contact points entered yet.")
        print("Fill in CONTACT_POINTS dict with (x, y) from ParaView at z=0.008.")

    # Create plot
    out_file = os.path.join(OUTPUT_DIR, "angle_vs_radii.png")
    create_plot(sim_nr, sim_ang, out_file)
