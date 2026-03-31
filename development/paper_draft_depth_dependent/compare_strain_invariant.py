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
# CASE 1F: Constant Material Properties
# ==================================================
def compute_case1f():
    # Physical constants
    density_fluid = 1000  # kg/m^3
    rho = 2670  # kg/m^3
    g = 9.8  # m/s^2

    xi_o = -1.0
    d_0 = 0
    gamma_r = 35e9
    # Elastic constants
    mu = 32.04e9 + xi_o * d_0 * gamma_r  # shear modulus (Pa)
    lmbda = 32.04e9  # Lame's first parameter (Pa)

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

    # Piecewise definitions with tapering applied to deviatoric stress
    sigma_xx = Omega * (b_xx * (sigma_zz + Pf) - Pf) + (1 - Omega) * sigma_zz
    sigma_yy = Omega * (b_yy * (sigma_zz + Pf) - Pf) + (1 - Omega) * sigma_zz
    sigma_xy = Omega * (b_xy * (sigma_zz + Pf))

    # Convert total stress to effective stress
    sigma_xx = sigma_xx + Pf
    sigma_yy = sigma_yy + Pf
    sigma_zz = sigma_zz + Pf

    # Assemble stress tensor components
    stress = np.zeros((len(depths), 3, 3))
    stress[:, 0, 0] = sigma_xx
    stress[:, 1, 1] = sigma_yy
    stress[:, 2, 2] = sigma_zz
    stress[:, 0, 1] = sigma_xy
    stress[:, 1, 0] = sigma_xy

    # Inverse Hooke's law: epsilon_ij = 1/(2μ) s_ij - λ/(2μ(3λ+2μ)) s_kk δ_ij
    trace_coeff = -lmbda / (2 * mu * (3 * lmbda + 2 * mu))
    strain = np.zeros_like(stress)
    enum = np.trace(stress, axis1=1, axis2=2)  # vector of s_kk
    for i in range(3):
        for j in range(3):
            strain[:, i, j] = stress[:, i, j] / (2 * mu)
    for k in range(3):
        strain[:, k, k] += trace_coeff * enum

    # Compute invariants: I1, I2, xi
    I1 = np.trace(strain, axis1=1, axis2=2)
    I2 = np.einsum("nij,nij->n", strain, strain)
    xi = I1 / np.sqrt(I2)

    return depths, xi


# ==================================================
# CASE 3F: Depth-Varying Material Properties (TPV32)
# ==================================================
def compute_case3f():
    # Physical constants
    density_fluid = 1000  # kg/m^3
    g = 9.8  # m/s^2

    # Elastic constants (depth-varying via TPV32 profile)
    depth_knots = np.array(
        [0, 500, 1000, 1600, 2400, 3600, 5000, 9000, 11000, 15000], dtype=float
    )
    vp_knots = np.array(
        [2200, 3000, 3600, 4400, 4800, 5250, 5500, 5750, 6100, 6300], dtype=float
    )
    vs_knots = np.array(
        [1050, 1400, 1950, 2500, 2800, 3100, 3250, 3450, 3600, 3700], dtype=float
    )
    rho_knots = np.array(
        [2200, 2450, 2550, 2600, 2600, 2620, 2650, 2720, 2750, 2900], dtype=float
    )

    # Depth array (m)
    depths = np.linspace(0, 20000, 600)  # from surface to 20 km

    # Interpolate seismic properties
    vp = np.interp(depths, depth_knots, vp_knots, left=vp_knots[0], right=vp_knots[-1])
    vs = np.interp(depths, depth_knots, vs_knots, left=vs_knots[0], right=vs_knots[-1])
    rho_depth = np.interp(
        depths, depth_knots, rho_knots, left=rho_knots[0], right=rho_knots[-1]
    )

    # Derived elastic constants (arrays)
    mu_depth = rho_depth * vs**2
    lmbda_depth = np.maximum(rho_depth * vp**2 - 2.0 * mu_depth, 0.0)

    # Pore pressure
    Pf = density_fluid * g * depths

    # Overburden from depth-varying rock density
    rho_avg = 0.5 * (rho_depth[1:] + rho_depth[:-1])
    dz = np.diff(depths)
    overburden = np.concatenate(([0.0], np.cumsum(rho_avg * dz)))
    sigma_zz = -g * overburden  # compressive negative

    # Coefficients for horizontal and shear stresses
    b_xx = 0.926793
    b_yy = 1.073206
    b_xy = -0.8

    # Tapering coefficient Omega(depth)
    Omega = np.ones_like(depths)
    mask_taper = (depths > 15000) & (depths <= 20000)
    Omega[mask_taper] = (20000 - depths[mask_taper]) / 5000
    Omega[depths > 20000] = 0.0

    # Deviatoric stress with tapering
    sigma_xx = Omega * (b_xx * (sigma_zz + Pf) - Pf) + (1 - Omega) * sigma_zz
    sigma_yy = Omega * (b_yy * (sigma_zz + Pf) - Pf) + (1 - Omega) * sigma_zz
    sigma_xy = Omega * (b_xy * (sigma_zz + Pf))

    # Convert total stress to effective stress
    sigma_xx = sigma_xx + Pf
    sigma_yy = sigma_yy + Pf
    sigma_zz = sigma_zz + Pf

    # Assemble stress tensor components
    stress = np.zeros((len(depths), 3, 3))
    stress[:, 0, 0] = sigma_xx
    stress[:, 1, 1] = sigma_yy
    stress[:, 2, 2] = sigma_zz
    stress[:, 0, 1] = sigma_xy
    stress[:, 1, 0] = sigma_xy

    # Inverse Hooke's law with depth-varying μ, λ
    trace_coeff = -lmbda_depth / (2 * mu_depth * (3 * lmbda_depth + 2 * mu_depth))
    strain = stress / (2.0 * mu_depth[:, None, None])  # first term
    enum = np.trace(stress, axis1=1, axis2=2)  # s_kk
    for k in range(3):
        strain[:, k, k] += trace_coeff * enum  # add diagonal second term

    # Invariants: I1, I2, xi
    I1 = np.trace(strain, axis1=1, axis2=2)
    I2 = np.einsum("nij,nij->n", strain, strain)
    I2 = np.maximum(I2, 1e-18)  # stabilize
    xi = I1 / np.sqrt(I2)

    return depths, xi


# ==================================================
# Main: Compute and Plot Comparison
# ==================================================
if __name__ == "__main__":
    # Compute strain invariant ratio for both cases
    depths_1f, xi_1f = compute_case1f()
    depths_3f, xi_3f = compute_case3f()

    # Create comparison plot
    fig, ax = plt.subplots(figsize=(8, 10))
    tick_prop = mpl.font_manager.FontProperties(family="DejaVu Sans", size=16)

    # Plot both cases
    ax.plot(
        xi_1f,
        depths_1f / 1e3,
        label="Case 1f (Constant Properties)",
        linewidth=4,
        color=THEME_ACCENTS[1],
    )
    ax.plot(
        xi_3f,
        depths_3f / 1e3,
        label="Case 3f (Depth-Varying TPV32)",
        linewidth=4,
        color=THEME_ACCENTS[2],
        linestyle="--",
    )

    # Format plot
    ax.invert_yaxis()
    ax.set_ylabel("Depth (km)", fontsize=22)
    ax.set_xlabel(r"Strain Invariant Ratio $\xi = I_1/\sqrt{I_2}$", fontsize=22)
    ax.set_title(
        "Comparison of Strain Invariant Ratio:\nCase 1f vs Case 3f", fontsize=24
    )
    # ax.legend(loc="best", fontsize=16, framealpha=0.9)
    # ax.grid(True, which="both", linestyle=":", alpha=0.6)

    # Set tick properties
    for lab in ax.get_xticklabels() + ax.get_yticklabels():
        lab.set_fontproperties(tick_prop)

    # Save figure
    plt.tight_layout()
    plt.savefig("strain_invariant_comparison_case1f_vs_case3f.png", dpi=300)
    print("Plot saved as: strain_invariant_comparison_case1f_vs_case3f.png")
    plt.show()

    # Print some statistics
    print("\n" + "=" * 60)
    print("Strain Invariant Ratio (ξ) Statistics")
    print("=" * 60)
    print(f"Case 1f - Min: {np.min(xi_1f):.4f}, Max: {np.max(xi_1f):.4f}")
    print(f"Case 3f - Min: {np.min(xi_3f):.4f}, Max: {np.max(xi_3f):.4f}")
    print(f"Maximum difference: {np.max(np.abs(xi_1f - xi_3f)):.4f}")
    print("=" * 60)
