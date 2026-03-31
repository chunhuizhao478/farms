"""Compare stress drop between case1f (hydrostatic) and case2f (overpressure).

Stress drop is defined as:
    Δτ = |σ_xy| - μ_d * |σ_n|

where σ_n is the effective normal stress (σ_yy in effective stress form).
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler

# ----------------------------------------------
# THEME (consistent with case files)
# ----------------------------------------------
THEME_ACCENTS = ["#000000", "#577399", "#C22B48", "#FFAB40", "#579499", "#5C5799"]
THEME_LT1 = "#FFFFFF"
THEME_DK1 = "#000000"

mpl.rcParams.update(
    {
        "axes.prop_cycle": cycler("color", THEME_ACCENTS),
        "figure.facecolor": THEME_LT1,
        "axes.facecolor": THEME_LT1,
        "savefig.facecolor": THEME_LT1,
        "text.color": THEME_DK1,
        "axes.labelcolor": THEME_DK1,
        "xtick.color": THEME_DK1,
        "ytick.color": THEME_DK1,
        "axes.edgecolor": THEME_DK1,
        "grid.color": THEME_DK1,
        "legend.edgecolor": THEME_DK1,
    }
)


# ==================================================
# CASE 1F: Hydrostatic Pore Pressure
# ==================================================
def compute_case1f():
    """Case 1f: Constant material properties with hydrostatic pore pressure."""
    # Physical constants
    density_fluid = 1000  # kg/m^3
    rho = 2670  # kg/m^3
    g = 9.8  # m/s^2

    # Friction coefficient
    mu_d = 0.6  # dynamic friction

    # Depth array (m)
    depths = np.linspace(0, 20000, 600)  # from surface to 20 km

    # Pore pressure and vertical stress
    Pf = density_fluid * g * depths
    sigma_zz = -rho * g * depths

    # Coefficients for horizontal and shear stresses
    b_xx = 0.926793
    b_yy = 1.073206
    b_xy = -0.8

    # Tapering coefficient Omega(depth)
    Omega = np.ones_like(depths)
    Omega[(depths > 15000) & (depths <= 20000)] = (
        20000 - depths[(depths > 15000) & (depths <= 20000)]
    ) / 5000
    Omega[depths > 20000] = 0.0

    # Total stresses with tapering
    sigma_xx = Omega * (b_xx * (sigma_zz + Pf) - Pf) + (1 - Omega) * sigma_zz
    sigma_yy = Omega * (b_yy * (sigma_zz + Pf) - Pf) + (1 - Omega) * sigma_zz
    sigma_xy = Omega * (b_xy * (sigma_zz + Pf))

    # Convert to effective stress
    sigma_xx_eff = sigma_xx + Pf
    sigma_yy_eff = sigma_yy + Pf
    sigma_zz_eff = sigma_zz + Pf

    # Stress drop: Δτ = |σ_xy| - μ_d * |σ_n|
    # Note: σ_xy is already in effective form (no Pf term in shear)
    stress_drop = np.abs(sigma_xy) - mu_d * np.abs(sigma_yy_eff)

    return depths, stress_drop, sigma_xy, sigma_yy_eff, Pf


# ==================================================
# CASE 2F: Overpressure with Depth-Dependent Transition
# ==================================================
def compute_case2f():
    """Case 2f: Constant material properties with overpressure transition."""
    # Physical constants
    density_fluid = 1000  # kg/m^3
    rho = 2670  # kg/m^3
    g = 9.8  # m/s^2

    # Friction coefficient
    mu_d = 0.6  # dynamic friction

    # Transition zone parameters
    A = 6000.0  # depth where transition begins (m)
    B = 8000.0  # depth where transition ends (m)

    # Depth array (m)
    depths = np.linspace(0, 20000, 600)  # from surface to 20 km

    # Pore pressure (piecewise with transition to overpressure)
    delta_rho = rho - density_fluid
    Pf = np.empty_like(depths)

    # Region 1: Hydrostatic above A
    mask1 = depths <= A
    Pf[mask1] = density_fluid * g * depths[mask1]

    # Region 2: Linear-gradient transition A < z ≤ B
    mask2 = (depths > A) & (depths <= B)
    z2 = depths[mask2]
    Pf_A = density_fluid * g * A
    Pf[mask2] = Pf_A + g * (
        density_fluid * (z2 - A) + 0.5 * delta_rho * (z2 - A) ** 2 / (B - A)
    )

    # Region 3: Over-pressured below B
    mask3 = depths > B
    z3 = depths[mask3]
    Pf_B = density_fluid * g * B + 0.5 * g * delta_rho * (B - A)
    Pf[mask3] = Pf_B + rho * g * (z3 - B)

    # Vertical stress
    sigma_zz = -rho * g * depths

    # Coefficients for horizontal and shear stresses
    b_xx = 0.926793
    b_yy = 1.073206
    b_xy = -0.8

    # Tapering coefficient Omega(depth)
    Omega = np.ones_like(depths)
    Omega[(depths > 15000) & (depths <= 20000)] = (
        20000 - depths[(depths > 15000) & (depths <= 20000)]
    ) / 5000
    Omega[depths > 20000] = 0.0

    # Total stresses with tapering
    sigma_xx = Omega * (b_xx * (sigma_zz + Pf) - Pf) + (1 - Omega) * sigma_zz
    sigma_yy = Omega * (b_yy * (sigma_zz + Pf) - Pf) + (1 - Omega) * sigma_zz
    sigma_xy = Omega * (b_xy * (sigma_zz + Pf))

    # Convert to effective stress
    sigma_xx_eff = sigma_xx + Pf
    sigma_yy_eff = sigma_yy + Pf
    sigma_zz_eff = sigma_zz + Pf

    # Stress drop: Δτ = |σ_xy| - μ_d * |σ_n|
    stress_drop = np.abs(sigma_xy) - mu_d * np.abs(sigma_yy_eff)

    return depths, stress_drop, sigma_xy, sigma_yy_eff, Pf


# ==================================================
# Main: Compute and Plot Comparison
# ==================================================
if __name__ == "__main__":
    # Compute stress drop for both cases
    depths_1f, stress_drop_1f, sigma_xy_1f, sigma_yy_1f, Pf_1f = compute_case1f()
    depths_2f, stress_drop_2f, sigma_xy_2f, sigma_yy_2f, Pf_2f = compute_case2f()

    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 10))
    tick_prop = mpl.font_manager.FontProperties(family="DejaVu Sans", size=16)

    # ===== Left subplot: Stress Drop Comparison =====
    ax1.plot(
        stress_drop_1f / 1e6,
        depths_1f / 1e3,
        label="Case 1 (Hydrostatic)",
        linewidth=4,
        color=THEME_ACCENTS[1],
    )
    ax1.plot(
        stress_drop_2f / 1e6,
        depths_2f / 1e3,
        label="Case 2 (Overpressure)",
        linewidth=4,
        color=THEME_ACCENTS[2],
        linestyle="--",
    )

    # Add vertical line at zero
    ax1.axvline(x=0, color=THEME_DK1, linestyle=":", alpha=0.5, linewidth=1)

    # Format left subplot
    ax1.invert_yaxis()
    ax1.set_ylabel("Depth (km)", fontsize=22)
    ax1.set_xlabel(
        r"Stress Drop $\Delta\tau = |\sigma_{xy}| - \mu_d |\sigma_n|$ (MPa)",
        fontsize=22,
    )
    ax1.set_title("Stress Drop Comparison", fontsize=24)
    ax1.legend(loc="best", fontsize=18, framealpha=0.9)
    ax1.grid(True, which="both", linestyle=":", alpha=0.6)
    for lab in ax1.get_xticklabels() + ax1.get_yticklabels():
        lab.set_fontproperties(tick_prop)

    # ===== Right subplot: Pore Pressure Comparison =====
    ax2.plot(
        Pf_1f / 1e6,
        depths_1f / 1e3,
        label="Case 1 (Hydrostatic)",
        linewidth=4,
        color=THEME_ACCENTS[1],
    )
    ax2.plot(
        Pf_2f / 1e6,
        depths_2f / 1e3,
        label="Case 2 (Overpressure)",
        linewidth=4,
        color=THEME_ACCENTS[2],
        linestyle="--",
    )

    # Highlight transition zone for case 2f
    ax2.axhspan(6, 8, alpha=0.2, color=THEME_ACCENTS[2], label="Transition Zone")

    # Format right subplot
    ax2.invert_yaxis()
    ax2.set_ylabel("Depth (km)", fontsize=22)
    ax2.set_xlabel(r"Pore Pressure $P_f$ (MPa)", fontsize=22)
    ax2.set_title("Pore Pressure Profiles", fontsize=24)
    ax2.legend(loc="best", fontsize=18, framealpha=0.9)
    ax2.grid(True, which="both", linestyle=":", alpha=0.6)
    for lab in ax2.get_xticklabels() + ax2.get_yticklabels():
        lab.set_fontproperties(tick_prop)

    # Save figure
    plt.tight_layout()
    plt.savefig("stress_drop_comparison_case1f_vs_case2f.png", dpi=300)
    print("Plot saved as: stress_drop_comparison_case1f_vs_case2f.png")
    plt.show()

    # ===== Print Statistics =====
    print("\n" + "=" * 70)
    print("Stress Drop (Δτ) Statistics")
    print("=" * 70)
    print(f"Case 1f (Hydrostatic):")
    print(
        f"  Min: {np.min(stress_drop_1f) / 1e6:8.2f} MPa at {depths_1f[np.argmin(stress_drop_1f)] / 1e3:.2f} km"
    )
    print(
        f"  Max: {np.max(stress_drop_1f) / 1e6:8.2f} MPa at {depths_1f[np.argmax(stress_drop_1f)] / 1e3:.2f} km"
    )
    print(f"\nCase 2f (Overpressure):")
    print(
        f"  Min: {np.min(stress_drop_2f) / 1e6:8.2f} MPa at {depths_2f[np.argmin(stress_drop_2f)] / 1e3:.2f} km"
    )
    print(
        f"  Max: {np.max(stress_drop_2f) / 1e6:8.2f} MPa at {depths_2f[np.argmax(stress_drop_2f)] / 1e3:.2f} km"
    )
    print(
        f"\nMaximum difference: {np.max(np.abs(stress_drop_1f - stress_drop_2f)) / 1e6:.2f} MPa"
    )
    print("=" * 70)

    # Additional analysis: Where stress drop becomes negative (locked fault)
    locked_1f = depths_1f[stress_drop_1f < 0]
    locked_2f = depths_2f[stress_drop_2f < 0]

    if len(locked_1f) > 0:
        print(
            f"\nCase 1f: Fault locked (Δτ < 0) from {locked_1f[0] / 1e3:.2f} km to {locked_1f[-1] / 1e3:.2f} km"
        )
    else:
        print(f"\nCase 1f: Fault not locked at any depth (Δτ ≥ 0 everywhere)")

    if len(locked_2f) > 0:
        print(
            f"Case 2f: Fault locked (Δτ < 0) from {locked_2f[0] / 1e3:.2f} km to {locked_2f[-1] / 1e3:.2f} km"
        )
    else:
        print(f"Case 2f: Fault not locked at any depth (Δτ ≥ 0 everywhere)")
    print("=" * 70)
